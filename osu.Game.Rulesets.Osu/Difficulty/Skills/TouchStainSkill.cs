// Copyright (c) ppy Pty Ltd <contact@ppy.sh>. Licensed under the MIT Licence.
// See the LICENCE file in the repository root for full licence text.

using System;
using osu.Game.Rulesets.Difficulty.Preprocessing;
using osu.Game.Rulesets.Difficulty.Utils;
using osu.Game.Rulesets.Mods;
using osu.Game.Rulesets.Objects;
using osu.Game.Rulesets.Osu.Difficulty.Preprocessing;
using osu.Game.Rulesets.Osu.Objects;

namespace osu.Game.Rulesets.Osu.Difficulty.Skills
{
    /// <summary>
    /// Represents the skill required to correctly aim and tap to the notes in a beatmap, using two hands on a touchscreen.
    /// </summary>
    public class TouchStrainSkill : OsuStrainSkill
    {
        protected override int HistoryLength => 60;
        private const int sequence_length = 4;

        // The actions utilized are bitmasked into an integer
        private const int l = 0;
        private const int r = 1;

        private const double aim_multiplier = 23.25;
        private const double aim_decay = 0.15;

        private const double speed_multiplier = 1375;
        private const double speed_decay = 0.3;

        private const double singletap_adjust = 1.7; // Difficulty multiplier between singletapping and alternating

        private const double overall_multiplier = 0.97; // Multiplier for final value

        // Constants for hand-coordination bonuses
        private const double coordination_aim_max_bonus = 0.62;
        private const double coordination_speed_max_bonus = 0.18; // Exact constant such that fully alternating between two hands rewards as much as streams

        // We keep an array of actions for convenience in code
        private readonly int[] actions =
        {
            l, r
        };

        private readonly double clockRate;
        private readonly double greatWindow;
        private readonly bool withSliders;
        private readonly bool calculatingAim;

        // Having the weighted strain of the most recent note stored globally to determine the initial strain of each strain section
        private double weightedAimStrain = 1;
        private double weightedSpeedStrain = 1;
        private double weightedRhythm = 1;

        // The aim, speed, and likelihood of each sequence are stored in an array with each index representing the bitmask of the sequence
        private readonly double[] sequenceAimStrain = new double[1 << sequence_length];
        private readonly double[] sequenceSpeedStrain = new double[1 << sequence_length];
        private readonly double[] likelihood = new double[1 << sequence_length]; // Should add up to 1.0
        private readonly long[] fullSequence = new long[1 << sequence_length]; // Contains full sequence up to 63 notes to maintain approximate precision

        // Constant table for fast base-2 logarithm calculation
        private readonly int[] logTable =
        {
            63, 0, 58, 1, 59, 47, 53, 2,
            60, 39, 48, 27, 54, 33, 42, 3,
            61, 51, 37, 40, 49, 18, 28, 20,
            55, 30, 34, 11, 43, 14, 22, 4,
            62, 57, 46, 52, 38, 26, 32, 41,
            50, 36, 17, 19, 29, 10, 13, 21,
            56, 45, 25, 31, 35, 16, 9, 12,
            44, 24, 15, 8, 23, 7, 6, 5
        };

        public TouchStrainSkill(Mod[] mods, double clockRate, double hitWindowGreat, bool withSliders, bool calculatingAim)
            : base(mods)
        {
            greatWindow = hitWindowGreat;
            this.withSliders = withSliders;
            this.clockRate = clockRate;
            this.calculatingAim = calculatingAim;

            // Setting the initial strain of each sequence to 1
            for (int i = 0; i < 1 << sequence_length; i++)
            {
                sequenceAimStrain[i] = 1;
                sequenceSpeedStrain[i] = 1;
                fullSequence[i] = i;
            }

            // The likelihood of hitting the first note with either hand is 50%
            likelihood[l] = 0.5;
            likelihood[r] = 0.5;
        }

        protected override double CalculateInitialStrain(double time)
        {
            double timeDifference = time - Previous[0].StartTime;
            if (calculatingAim)
                return overall_multiplier * weightedAimStrain * decay(aim_decay, timeDifference);

            return overall_multiplier * weightedRhythm * weightedSpeedStrain * decay(speed_decay, timeDifference);
        }

        protected override double StrainValueAt(DifficultyHitObject current)
        {
            weightedAimStrain = 0;
            weightedSpeedStrain = 0;
            weightedRhythm = 0;

            double[] newAimStrain = new double[1 << sequence_length];
            double[] newSpeedStrain = new double[1 << sequence_length];
            double[] newlikelihood = new double[1 << sequence_length];
            long[] newFullSequence = new long[1 << sequence_length];
            int currentSequenceLength = Math.Min(sequence_length, Previous.Count + 1);
            int currentHistoryLength = Math.Min(HistoryLength, Previous.Count + 1);

            HitObject[] trueHistory = new HitObject[currentHistoryLength];
            // Raw HitObjects as opposed to DifficultyHitObjects will be required to reflect true notes hit
            trueHistory[0] = current.LastObject;

            for (int j = 1; j < currentHistoryLength; j++)
            {
                var previous = Previous[j - 1].LastObject;
                trueHistory[j] = previous;
            }

            // Performing a weighted sum of all sequences, simulating hitting the current note with either the left or right hand
            for (int i = 0; i < 1 << currentSequenceLength; i++)
            {
                // The strain of the sequence being added to is decayed
                double prevAim = sequenceAimStrain[i] * decay(aim_decay, current.DeltaTime);
                double prevSpeed = sequenceSpeedStrain[i] * decay(speed_decay, current.DeltaTime);
                // The raw aim and speed at the current note for each action is stored
                double[] rawAim = new double[2];
                double[] rawSpeed = new double[2];
                double[] rawRhythm = new double[2];

                long approxExtension = fullSequence[i];

                // What happens when we hit the current note with either hand, given a preexisting "history" of hand choices
                foreach (int hand in actions)
                {
                    long bitFixedExtension = approxExtension; // A version of approxExtension such that all bits that are the same as hand are set to 1
                    long invert = (1L << currentHistoryLength) - 1;
                    if (hand == 0)
                        bitFixedExtension ^= invert;

                    // Identifying relevant HitObject indexes from the bitmask sequence
                    HitObject[] lastSame = new HitObject[HistoryLength];
                    HitObject lastSwap = null;

                    int lastSwapIndex = rightmostBitLocation(bitFixedExtension ^ invert);
                    if (lastSwapIndex < currentHistoryLength)
                        lastSwap = trueHistory[lastSwapIndex];

                    int mostRecentSame = 0;

                    while (mostRecentSame < HistoryLength)
                    {
                        long leastSignificantBit = (bitFixedExtension & (-bitFixedExtension));
                        int lastSameIndex = rightmostBitLocation(bitFixedExtension);

                        if (lastSameIndex >= currentHistoryLength) break;

                        lastSame[mostRecentSame] = trueHistory[lastSameIndex];
                        mostRecentSame++;
                        bitFixedExtension -= leastSignificantBit;
                    }

                    if (mostRecentSame > 0)
                    {
                        // Calculating the individual strains
                        DifficultyHitObject simulatedCurrent = new OsuDifficultyHitObject(current.BaseObject, lastSame[1], lastSame[0], clockRate);
                        ReverseQueue<DifficultyHitObject> simulatedPrevious = new ReverseQueue<DifficultyHitObject>(HistoryLength);

                        const int min_required = 3; // Minimum DifficultyHitObjects required for aim/speed calculation, the rest will be "lazily" generated, only consisting of deltaTimes

                        for (int j = mostRecentSame - 2; j > min_required; j--)
                        {
                            simulatedPrevious.Enqueue(new OsuDifficultyHitObject(lastSame[j], j + 2 == mostRecentSame ? null : lastSame[j + 2], lastSame[j + 1], clockRate, false));
                        }

                        for (int j = Math.Min(min_required, mostRecentSame - 2); j >= 0; j--)
                        {
                            simulatedPrevious.Enqueue(new OsuDifficultyHitObject(lastSame[j], lastSame[j + 2], lastSame[j + 1], clockRate));
                        }

                        // Adding a hand-coordination bonus based on time since the last hand "swap"
                        double aimBonus = 1;
                        double speedBonus = 1;
                        double angleBonus = 1;

                        if (lastSwap != null)
                        {
                            // Increases aim coordination bonus if the most recent instance of the "other hand" is in between the current object and the previous object with the actual hand
                            var simulatedSwap = new OsuDifficultyHitObject(current.BaseObject, lastSame[0], lastSwap, clockRate);
                            double? angle = simulatedSwap.Angle;
                            if (angle != null)
                                angleBonus += 1 / (1 + Math.Pow(Math.E, -(angle.Value * 180 / Math.PI - 108) / 9));

                            double coordinationFactor = simulatedCurrent.DeltaTime / (simulatedCurrent.DeltaTime + simulatedSwap.DeltaTime);
                            coordinationFactor = Math.Min(1, 3.5 * Math.Pow(coordinationFactor, 3));
                            aimBonus += coordination_aim_max_bonus * coordinationFactor;
                            speedBonus += coordination_speed_max_bonus * coordinationFactor;
                        }

                        rawAim[hand] = aim_multiplier * Math.Pow(aimBonus * calcAim(simulatedCurrent, simulatedPrevious), angleBonus);
                        rawSpeed[hand] = speed_multiplier * speedBonus * calcSpeed(simulatedCurrent, simulatedPrevious);
                        rawRhythm[hand] = calculateRhythmBonus(simulatedCurrent, simulatedPrevious);
                    }
                }

                double leftStrain = curveStrain(rawAim[l], rawRhythm[l] * rawSpeed[l]);
                double rightStrain = curveStrain(rawAim[r], rawRhythm[r] * rawSpeed[r]);

                // Chance of hitting the current note with either action is a ratio
                double[] chances = new double[2];
                chances[l] = Math.Abs(leftStrain + rightStrain) <= double.Epsilon ? 0.5 : curveProbability(rightStrain / (leftStrain + rightStrain));
                chances[r] = 1.0 - chances[l];

                foreach (int hand in actions)
                {
                    int newSequence = appendAction(i, hand);

                    double chance = chances[hand] * likelihood[i];

                    if (chance > newlikelihood[newSequence]) // Set the full sequence to the more likely of the two cut sequences that "merge" into the newSequence
                        newFullSequence[newSequence] = appendHistory(fullSequence[i], hand);

                    newlikelihood[newSequence] += chance;

                    // Pushing the strain to the new sequence, starting from the current note
                    newAimStrain[newSequence] += chance * (prevAim + rawAim[hand]);
                    newSpeedStrain[newSequence] += chance * (prevSpeed + rawSpeed[hand]);
                    // Pushing the strain to the overall weighted aim and speed of the current note
                    weightedAimStrain += chance * (prevAim + rawAim[hand]);
                    weightedSpeedStrain += chance * (prevSpeed + rawSpeed[hand]);
                    weightedRhythm += chance * rawRhythm[hand];
                }
            }

            int newSequenceLength = Math.Min(sequence_length, currentSequenceLength + 1);

            // Dividing by the likelihood to return the weighted average of all past sequences that contributed to the current one
            for (int i = 0; i < 1 << newSequenceLength; i++)
            {
                if (newlikelihood[i] > double.Epsilon)
                {
                    newAimStrain[i] /= newlikelihood[i];
                    newSpeedStrain[i] /= newlikelihood[i];
                }
            }

            // Updating the cached arrays
            Array.Copy(newAimStrain, sequenceAimStrain, 1 << sequence_length);
            Array.Copy(newSpeedStrain, sequenceSpeedStrain, 1 << sequence_length);
            Array.Copy(newlikelihood, likelihood, 1 << sequence_length);
            Array.Copy(newFullSequence, fullSequence, 1 << sequence_length);

            // Curving the aim and speed strains to better match old values
            if (calculatingAim)
                return overall_multiplier * weightedAimStrain;

            return overall_multiplier * weightedRhythm * weightedSpeedStrain;
        }

        private double curveProbability(double chance) => chance <= 0.5 ? Math.Pow(2 * chance, 4) / 2 : 1 - Math.Pow(2 * (1 - chance), 4) / 2;

        private double decay(double decayBase, double ms) => Math.Pow(decayBase, ms / 1000);

        private double curveStrain(double aimStrain, double speedStrain) => Math.Pow(Math.Pow(aimStrain, 3.0 / 2.0) + Math.Pow(speedStrain, 3.0 / 2.0), 2.0 / 3.0);

        private int appendAction(int oldSequence, int toAppend) => ((oldSequence << 1) | toAppend) & ((1 << sequence_length) - 1);

        private long appendHistory(long oldSequence, long toAppend) => ((oldSequence << 1) | toAppend) & ((1 << HistoryLength) - 1);

        private int rightmostBitLocation(long value)
        {
            return fastLog2((ulong)(value & -value));
        }

        private int fastLog2(ulong value) => logTable[(int)((value * 0x07EDD5E59A4E28C2) >> 58)];

        /// <summary>
        /// Calculates a rhythm multiplier for the difficulty of the tap associated with historic data of the current <see cref="OsuDifficultyHitObject"/>.
        /// </summary>
        private double calculateRhythmBonus(DifficultyHitObject current, ReverseQueue<DifficultyHitObject> previous)
        {
            const double rhythm_multiplier = 0.75;
            const int history_time_max = 5000; // 5 seconds of calculatingRhythmBonus max.
            if (current.BaseObject is Spinner)
                return 0;

            int previousIslandSize = 0;

            double rhythmComplexitySum = 0;
            int islandSize = 1;
            double startRatio = 0; // store the ratio of the current start of an island to buff for tighter rhythms

            bool firstDeltaSwitch = false;

            int rhythmStart = 0;

            while (rhythmStart < previous.Count - 2 && current.StartTime - previous[rhythmStart].StartTime < history_time_max)
                rhythmStart++;

            for (int i = rhythmStart; i > 0; i--)
            {
                OsuDifficultyHitObject currObj = (OsuDifficultyHitObject)previous[i - 1];
                OsuDifficultyHitObject prevObj = (OsuDifficultyHitObject)previous[i];
                OsuDifficultyHitObject lastObj = (OsuDifficultyHitObject)previous[i + 1];

                double currHistoricalDecay = (history_time_max - (current.StartTime - currObj.StartTime)) / history_time_max; // scales note 0 to 1 from history to now

                currHistoricalDecay = Math.Min((double)(previous.Count - i) / previous.Count, currHistoricalDecay); // either we're limited by time or limited by object count.

                double currDelta = currObj.StrainTime;
                double prevDelta = prevObj.StrainTime;
                double lastDelta = lastObj.StrainTime;
                double currRatio = 1.0 + 6.0 * Math.Min(0.5,
                    Math.Pow(Math.Sin(Math.PI / (Math.Min(prevDelta, currDelta) / Math.Max(prevDelta, currDelta))), 2)); // fancy function to calculate rhythmbonuses.

                double windowPenalty = Math.Min(1, Math.Max(0, Math.Abs(prevDelta - currDelta) - greatWindow * 0.6) / (greatWindow * 0.6));

                windowPenalty = Math.Min(1, windowPenalty);

                double effectiveRatio = windowPenalty * currRatio;

                if (firstDeltaSwitch)
                {
                    if (!(prevDelta > 1.25 * currDelta || prevDelta * 1.25 < currDelta))
                    {
                        if (islandSize < 7)
                            islandSize++; // island is still progressing, count size.
                    }
                    else
                    {
                        if (previous[i - 1].BaseObject is Slider) // bpm change is into slider, this is easy acc window
                            effectiveRatio *= 0.125;

                        if (previous[i].BaseObject is Slider) // bpm change was from a slider, this is easier typically than circle -> circle
                            effectiveRatio *= 0.25;

                        if (previousIslandSize == islandSize) // repeated island size (ex: triplet -> triplet)
                            effectiveRatio *= 0.25;

                        if (previousIslandSize % 2 == islandSize % 2) // repeated island polartiy (2 -> 4, 3 -> 5)
                            effectiveRatio *= 0.50;

                        if (lastDelta > prevDelta + 10 && prevDelta > currDelta + 10) // previous increase happened a note ago, 1/1->1/2-1/4, dont want to buff this.
                            effectiveRatio *= 0.125;

                        rhythmComplexitySum += Math.Sqrt(effectiveRatio * startRatio) * currHistoricalDecay * Math.Sqrt(4 + islandSize) / 2 * Math.Sqrt(4 + previousIslandSize) / 2;

                        startRatio = effectiveRatio;

                        previousIslandSize = islandSize; // log the last island size.

                        if (prevDelta * 1.25 < currDelta) // we're slowing down, stop counting
                            firstDeltaSwitch = false; // if we're speeding up, this stays true and  we keep counting island size.

                        islandSize = 1;
                    }
                }
                else if (prevDelta > 1.25 * currDelta) // we want to be speeding up.
                {
                    // Begin counting island until we change speed again.
                    firstDeltaSwitch = true;
                    startRatio = effectiveRatio;
                    islandSize = 1;
                }
            }

            return Math.Sqrt(4 + rhythmComplexitySum * rhythm_multiplier) / 2; //produces multiplier that can be applied to strain. range [1, infinity) (not really though)
        }

        /// <summary>
        /// Represents the skill required to correctly aim at every object in the map with one hand, using a uniform CircleSize and normalized distances.
        /// </summary>
        private double calcAim(DifficultyHitObject current, ReverseQueue<DifficultyHitObject> previous)
        {
            const double wide_angle_multiplier = 1.5;
            const double acute_angle_multiplier = 2.0;
            const double slider_multiplier = 1.5;
            const double velocity_change_multiplier = 0.75;

            if (current.BaseObject is Spinner || previous.Count <= 1 || previous[0].BaseObject is Spinner)
                return 0;

            var osuCurrObj = (OsuDifficultyHitObject)current;
            var osuLastObj = (OsuDifficultyHitObject)previous[0];
            var osuLastLastObj = (OsuDifficultyHitObject)previous[1];

            // Calculate the velocity to the current hitobject, which starts with a base distance / time assuming the last object is a hitcircle.
            double currVelocity = osuCurrObj.LazyJumpDistance / osuCurrObj.StrainTime;

            // But if the last object is a slider, then we extend the travel velocity through the slider into the current object.
            if (osuLastObj.BaseObject is Slider && withSliders)
            {
                double travelVelocity = osuLastObj.TravelDistance / osuLastObj.TravelTime; // calculate the slider velocity from slider head to slider end.
                double movementVelocity = osuCurrObj.MinimumJumpDistance / osuCurrObj.MinimumJumpTime; // calculate the movement velocity from slider end to current object

                currVelocity = Math.Max(currVelocity, movementVelocity + travelVelocity); // take the larger total combined velocity.
            }

            // As above, do the same for the previous hitobject.
            double prevVelocity = osuLastObj.LazyJumpDistance / osuLastObj.StrainTime;

            if (osuLastLastObj.BaseObject is Slider && withSliders)
            {
                double travelVelocity = osuLastLastObj.TravelDistance / osuLastLastObj.TravelTime;
                double movementVelocity = osuLastObj.MinimumJumpDistance / osuLastObj.MinimumJumpTime;

                prevVelocity = Math.Max(prevVelocity, movementVelocity + travelVelocity);
            }

            double wideAngleBonus = 0;
            double acuteAngleBonus = 0;
            double sliderBonus = 0;
            double velocityChangeBonus = 0;

            double aimStrain = currVelocity; // Start strain with regular velocity.

            if (Math.Max(osuCurrObj.StrainTime, osuLastObj.StrainTime) < 1.25 * Math.Min(osuCurrObj.StrainTime, osuLastObj.StrainTime)) // If rhythms are the same.
            {
                if (osuCurrObj.Angle != null && osuLastObj.Angle != null && osuLastLastObj.Angle != null)
                {
                    double currAngle = osuCurrObj.Angle.Value;
                    double lastAngle = osuLastObj.Angle.Value;
                    double lastLastAngle = osuLastLastObj.Angle.Value;

                    // Rewarding angles, take the smaller velocity as base.
                    double angleBonus = Math.Min(currVelocity, prevVelocity);

                    wideAngleBonus = calcWideAngleBonus(currAngle);
                    acuteAngleBonus = calcAcuteAngleBonus(currAngle);

                    if (osuCurrObj.StrainTime > 100) // Only buff deltaTime exceeding 300 bpm 1/2.
                        acuteAngleBonus = 0;
                    else
                    {
                        acuteAngleBonus *= calcAcuteAngleBonus(lastAngle) // Multiply by previous angle, we don't want to buff unless this is a wiggle type pattern.
                                           * Math.Min(angleBonus, 125 / osuCurrObj.StrainTime) // The maximum velocity we buff is equal to 125 / strainTime
                                           * Math.Pow(Math.Sin(Math.PI / 2 * Math.Min(1, (100 - osuCurrObj.StrainTime) / 25)), 2) // scale buff from 150 bpm 1/4 to 200 bpm 1/4
                                           * Math.Pow(Math.Sin(Math.PI / 2 * (Math.Clamp(osuCurrObj.LazyJumpDistance, 50, 100) - 50) / 50),
                                               2); // Buff distance exceeding 50 (radius) up to 100 (diameter).
                    }

                    // Penalize wide angles if they're repeated, reducing the penalty as the lastAngle gets more acute.
                    wideAngleBonus *= angleBonus * (1 - Math.Min(wideAngleBonus, Math.Pow(calcWideAngleBonus(lastAngle), 3)));
                    // Penalize acute angles if they're repeated, reducing the penalty as the lastLastAngle gets more obtuse.
                    acuteAngleBonus *= 0.5 + 0.5 * (1 - Math.Min(acuteAngleBonus, Math.Pow(calcAcuteAngleBonus(lastLastAngle), 3)));
                }
            }

            if (Math.Max(prevVelocity, currVelocity) != 0)
            {
                // We want to use the average velocity over the whole object when awarding differences, not the individual jump and slider path velocities.
                prevVelocity = (osuLastObj.LazyJumpDistance + osuLastLastObj.TravelDistance) / osuLastObj.StrainTime;
                currVelocity = (osuCurrObj.LazyJumpDistance + osuLastObj.TravelDistance) / osuCurrObj.StrainTime;

                // Scale with ratio of difference compared to 0.5 * max dist.
                double distRatio = Math.Pow(Math.Sin(Math.PI / 2 * Math.Abs(prevVelocity - currVelocity) / Math.Max(prevVelocity, currVelocity)), 2);

                // Reward for % distance up to 125 / strainTime for overlaps where velocity is still changing.
                double overlapVelocityBuff = Math.Min(125 / Math.Min(osuCurrObj.StrainTime, osuLastObj.StrainTime), Math.Abs(prevVelocity - currVelocity));

                // Reward for % distance slowed down compared to previous, paying attention to not award overlap
                double nonOverlapVelocityBuff = Math.Abs(prevVelocity - currVelocity)
                                                // do not award overlap
                                                * Math.Pow(Math.Sin(Math.PI / 2 * Math.Min(1, Math.Min(osuCurrObj.LazyJumpDistance, osuLastObj.LazyJumpDistance) / 100)), 2);

                // Choose the largest bonus, multiplied by ratio.
                velocityChangeBonus = Math.Max(overlapVelocityBuff, nonOverlapVelocityBuff) * distRatio;

                // Penalize for rhythm changes.
                velocityChangeBonus *= Math.Pow(Math.Min(osuCurrObj.StrainTime, osuLastObj.StrainTime) / Math.Max(osuCurrObj.StrainTime, osuLastObj.StrainTime), 2);
            }

            if (osuLastObj.TravelTime != 0)
            {
                // Reward sliders based on velocity.
                sliderBonus = osuLastObj.TravelDistance / osuLastObj.TravelTime;
            }

            // Add in acute angle bonus or wide angle bonus + velocity change bonus, whichever is larger.
            aimStrain += Math.Max(acuteAngleBonus * acute_angle_multiplier, wideAngleBonus * wide_angle_multiplier + velocityChangeBonus * velocity_change_multiplier);

            // Add in additional slider velocity bonus.
            if (withSliders)
                aimStrain += sliderBonus * slider_multiplier;

            return aimStrain;
        }

        /// <summary>
        /// Represents the skill required to tap with one finger with regards to keeping up with the speed at which objects need to be hit.
        /// </summary>
        private double calcSpeed(DifficultyHitObject current, ReverseQueue<DifficultyHitObject> previous)
        {
            const double single_spacing_threshold = 125;

            const double min_speed_bonus = 75; // ~200BPM
            const double speed_balancing_factor = 40;

            if (current.BaseObject is Spinner)
                return 0;

            // derive strainTime for calculation
            var osuCurrObj = (OsuDifficultyHitObject)current;
            var osuPrevObj = previous.Count > 0 ? (OsuDifficultyHitObject)previous[0] : null;

            double strainTime = osuCurrObj.StrainTime / singletap_adjust;
            double greatWindowFull = greatWindow * 2;

            // Cap deltatime to the OD 300 hitwindow.
            // 0.93 is derived from making sure 260bpm OD8 streams aren't nerfed harshly, whilst 0.92 limits the effect of the cap.
            strainTime /= Math.Clamp((strainTime / greatWindowFull) / 0.93, 0.92, 1);

            // derive speedBonus for calculation
            double speedBonus = 1.0;

            if (strainTime < min_speed_bonus)
                speedBonus = 1 + 0.75 * Math.Pow((min_speed_bonus - strainTime) / speed_balancing_factor, 2);

            double travelDistance = osuPrevObj?.TravelDistance ?? 0;
            double distance = Math.Min(single_spacing_threshold, travelDistance + osuCurrObj.MinimumJumpDistance);

            return (speedBonus + speedBonus * Math.Pow(distance / single_spacing_threshold, 3.5)) / strainTime;
        }

        private double calcWideAngleBonus(double angle) => Math.Pow(Math.Sin(3.0 / 4 * (Math.Min(5.0 / 6 * Math.PI, Math.Max(Math.PI / 6, angle)) - Math.PI / 6)), 2);

        private double calcAcuteAngleBonus(double angle) => 1 - calcWideAngleBonus(angle);
    }
}
