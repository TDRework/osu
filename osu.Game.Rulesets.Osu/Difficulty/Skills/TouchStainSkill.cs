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
        protected override int HistoryLength => sequence_length + min_previous;

        private readonly ReverseQueue<DifficultyHitObject>[] history =
        {
            new ReverseQueue<DifficultyHitObject>(hand_history_length),
            new ReverseQueue<DifficultyHitObject>(hand_history_length)
        };

        private const int sequence_length = 5;
        private const int hand_history_length = 32;
        private const int min_previous = 3;

        // The actions utilized are bitmasked into an integer
        private const int l = 0;
        private const int r = 1;

        private const double aim_multiplier = 26.25;
        private const double aim_decay = 0.15;

        private const double speed_multiplier = 1375;
        private const double speed_decay = 0.3;

        private const double singletap_adjust = 1.7; // Difficulty multiplier between singletapping and alternating

        private const double overall_multiplier = 0.92; // Multiplier for final value

        // Constants for hand-coordination bonuses
        private const double coordination_bonus_start = 400;
        private const double coordination_bonus_end = 75; // ~400BPM jump-speed
        private const double coordination_aim_max_bonus = 0.8;
        private const double coordination_speed_max_bonus = 0.2;
        private const double coordination_aim_decay_strength = 3;
        private const double coordination_speed_decay_strength = 3;

        // We keep an array of actions for convenience in code
        private readonly int[] actions =
        {
            l, r
        };

        private readonly double clockRate;
        private readonly double greatWindow;
        private readonly bool withSliders;

        // Having the weighted strain of the most recent note stored globally to determine the initial strain of each strain section
        private double weightedAimStrain = 1;
        private double weightedSpeedStrain = 1;
        private double currentRhythm;

        // The aim, speed, and likelihood of each sequence are stored in an array with each index representing the bitmask of the sequence
        private readonly double[] sequenceAimStrain = new double[1 << sequence_length];
        private readonly double[] sequenceSpeedStrain = new double[1 << sequence_length];
        private readonly double[] likelihood = new double[1 << sequence_length]; // Should add up to 1.0

        public TouchStrainSkill(Mod[] mods, double clockRate, double hitWindowGreat, bool withSliders)
            : base(mods)
        {
            greatWindow = hitWindowGreat;
            this.withSliders = withSliders;
            this.clockRate = clockRate;

            // Setting the initial strain of each sequence to 1
            for (int i = 0; i < 1 << sequence_length; i++)
            {
                sequenceAimStrain[i] = 1;
                sequenceSpeedStrain[i] = 1;
            }

            // The likelihood of hitting the first note with either hand is 50%
            likelihood[l] = 0.5;
            likelihood[r] = 0.5;
        }

        protected override double CalculateInitialStrain(double time)
        {
            double timeDifference = time - Previous[0].StartTime;
            return curveStrain(weightedAimStrain * decay(aim_decay, timeDifference),
                weightedSpeedStrain * currentRhythm * decay(speed_decay, timeDifference));
        }

        protected override double StrainValueAt(DifficultyHitObject current)
        {
            weightedAimStrain = 0;
            weightedSpeedStrain = 0;

            double[] newAimStrain = new double[1 << sequence_length];
            double[] newSpeedStrain = new double[1 << sequence_length];
            double[] newlikelihood = new double[1 << sequence_length];
            int currentSequenceLength = Math.Min(sequence_length, Previous.Count + 1);

            // Performing a weighted sum of all sequences, simulating hitting the current note with either the left or right hand
            for (int i = 0; i < 1 << currentSequenceLength; i++)
            {
                // The strain of the sequence being added to is decayed
                double prevAim = sequenceAimStrain[i] * decay(aim_decay, current.DeltaTime);
                double prevSpeed = sequenceSpeedStrain[i] * decay(speed_decay, current.DeltaTime);
                // The raw aim and speed at the current note for each action is stored
                double[] rawAim = new double[2];
                double[] rawSpeed = new double[2];

                foreach (var hand in actions)
                {
                    HitObject[] trueHistory = new HitObject[currentSequenceLength + min_previous];
                    trueHistory[0] = current.LastObject;

                    for (int j = 1; j < currentSequenceLength + min_previous; j++)
                    {
                        if (j <= Previous.Count)
                        {
                            var previous = Previous[j - 1].LastObject;
                            trueHistory[j] = previous;
                        }
                    }

                    // Identifying relevant HitObject indexes from the bitmask sequence
                    HitObject[] lastSame = new HitObject[min_previous];
                    HitObject lastSwap = null;
                    int mostRecentSame = 0;

                    for (int j = 0; j < currentSequenceLength; j++)
                    {
                        int lastHand = (i >> j) & 1;

                        if (mostRecentSame < min_previous && lastHand == hand)
                        {
                            lastSame[mostRecentSame] = trueHistory[j];
                            mostRecentSame++;
                        }
                        else if (lastSwap == null) lastSwap = trueHistory[j];
                    }

                    // Searching previous objects beyond the bitmask length if the requirement of three previous objects has not been met
                    int mostRecent = currentSequenceLength;

                    while (mostRecentSame < min_previous && mostRecent < Previous.Count)
                    {
                        lastSame[mostRecentSame] = trueHistory[mostRecent];
                        mostRecentSame++;
                        mostRecent++;
                    }

                    if (mostRecentSame != 0)
                    {
                        // Calculating the individual strains
                        DifficultyHitObject simulatedCurrent = new OsuDifficultyHitObject(current.BaseObject, lastSame[1], lastSame[0], clockRate);
                        DifficultyHitObject simulatedPrevious = lastSame[1] == null ? null : new OsuDifficultyHitObject(lastSame[0], lastSame[2], lastSame[1], clockRate);
                        DifficultyHitObject simulatedPreviousPrevious = lastSame[2] == null ? null : new OsuDifficultyHitObject(lastSame[1], null, lastSame[2], clockRate);
                        // Adding a hand-coordination bonus based on time since the last hand "swap"
                        double aimBonus = 1;
                        double speedBonus = 1;
                        double angleBonus = 1;

                        if (lastSwap != null)
                        {
                            // Increases aim coordination bonus if the most recent instance of the "other hand" is in between the current object and the previous object with the actual hand
                            var angle = new OsuDifficultyHitObject(current.BaseObject, lastSame[0], lastSwap, clockRate).Angle;

                            if (angle != null)
                            {
                                angleBonus += 1 / (1 + Math.Pow(Math.E, -(angle.Value * 180 / Math.PI - 108) / 9));
                            }

                            double ms = new OsuDifficultyHitObject(current.BaseObject, null, lastSwap, clockRate).StrainTime;
                            ms = Math.Min(Math.Max(ms, coordination_bonus_end), coordination_bonus_start);
                            aimBonus += coordination_aim_max_bonus * Math.Pow(coordination_bonus_start - ms, coordination_aim_decay_strength)
                                        / Math.Pow(coordination_bonus_start - coordination_bonus_end, coordination_aim_decay_strength);
                            speedBonus += coordination_speed_max_bonus * Math.Pow(coordination_bonus_start - ms, coordination_speed_decay_strength)
                                          / Math.Pow(coordination_bonus_start - coordination_bonus_end, coordination_speed_decay_strength);
                        }

                        rawAim[hand] = aimBonus * aim_multiplier * Math.Pow(calcAim(simulatedCurrent, simulatedPrevious, simulatedPreviousPrevious), angleBonus);
                        rawSpeed[hand] = speedBonus * speed_multiplier * calcSpeed(simulatedCurrent, simulatedPrevious);
                    }
                }

                double leftStrain = curveStrain(rawAim[l], rawSpeed[l]);
                double rightStrain = curveStrain(rawAim[r], rawSpeed[r]);

                // Chance of hitting the current note with either action is a ratio
                double[] chances = new double[2];
                chances[l] = Math.Abs(leftStrain + rightStrain) <= double.Epsilon ? 0.5 : curveProbability(rightStrain / (leftStrain + rightStrain));
                chances[r] = 1.0 - chances[l];

                foreach (var hand in actions)
                {
                    int newSequence = appendAction(i, hand);
                    double chance = chances[hand] * likelihood[i];
                    newlikelihood[newSequence] += chance;

                    // Pushing the strain to the new sequence, starting from the current note
                    newAimStrain[newSequence] += chance * (prevAim + rawAim[hand]);
                    newSpeedStrain[newSequence] += chance * (prevSpeed + rawSpeed[hand]);
                    // Pushing the strain to the overall weighted aim and speed of the current note
                    weightedAimStrain += chance * (prevAim + rawAim[hand]);
                    weightedSpeedStrain += chance * (prevSpeed + rawSpeed[hand]);
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

            double[] finalLikelihood = new double[2];

            // Identifing overall likelyhood of hitting the current note with either hand
            for (int i = 0; i < 1 << newSequenceLength; i++)
            {
                finalLikelihood[i & 1] += newlikelihood[i];
            }

            // Calculating rhythmic bonus based on the previous 32 notes most likely to be hit with the recent hand
            int recentHand = finalLikelihood[l] < finalLikelihood[r] ? r : l;
            var simulatedLast = history[recentHand].Count > 0 ? history[recentHand][0].BaseObject : null;
            var simulatedLastLast = history[recentHand].Count > 1 ? history[recentHand][1].BaseObject : null;
            var simulatedBase = simulatedLast == null ? current : new OsuDifficultyHitObject(current.BaseObject, simulatedLastLast, simulatedLast, clockRate);
            currentRhythm = calculateRhythmBonus(simulatedBase, recentHand);

            while (history[recentHand].Count > hand_history_length)
                history[recentHand].Dequeue();
            history[recentHand].Enqueue(simulatedBase);

            // Updating the cached arrays
            Array.Copy(newAimStrain, sequenceAimStrain, 1 << sequence_length);
            Array.Copy(newSpeedStrain, sequenceSpeedStrain, 1 << sequence_length);
            Array.Copy(newlikelihood, likelihood, 1 << sequence_length);

            // Curving the aim and speed strains to better match old values
            return overall_multiplier * curveStrain(weightedAimStrain, currentRhythm * weightedSpeedStrain);
        }

        private double curveProbability(double chance) => chance <= 0.5 ? Math.Pow(2 * chance, 4) / 2 : 1 - Math.Pow(2 * (1 - chance), 4) / 2;

        private double decay(double decayBase, double ms) => Math.Pow(decayBase, ms / 1000);

        private double curveStrain(double aimStrain, double speedStrain) => Math.Pow(Math.Pow(aimStrain, 3.0 / 2.0) + Math.Pow(speedStrain, 3.0 / 2.0), 2.0 / 3.0);

        private int appendAction(int oldSequence, int toAppend) => ((oldSequence << 1) | toAppend) & ((1 << sequence_length) - 1);

        /// <summary>
        /// Calculates a rhythm multiplier for the difficulty of the tap associated with historic data of the current <see cref="OsuDifficultyHitObject"/>.
        /// </summary>
        private double calculateRhythmBonus(DifficultyHitObject current, int hand)
        {
            const double rhythm_multiplier = 0.75;
            const int history_time_max = 5000; // 5 seconds of calculatingRhythmBonus max.
            if (current.BaseObject is Spinner)
                return 0;

            var previousHandNotes = history[hand];
            int previousIslandSize = 0;

            double rhythmComplexitySum = 0;
            int islandSize = 1;
            double startRatio = 0; // store the ratio of the current start of an island to buff for tighter rhythms

            bool firstDeltaSwitch = false;

            int rhythmStart = 0;

            while (rhythmStart < previousHandNotes.Count - 2 && current.StartTime - previousHandNotes[rhythmStart].StartTime < history_time_max)
                rhythmStart++;

            for (int i = rhythmStart; i > 0; i--)
            {
                OsuDifficultyHitObject currObj = (OsuDifficultyHitObject)previousHandNotes[i - 1];
                OsuDifficultyHitObject prevObj = (OsuDifficultyHitObject)previousHandNotes[i];
                OsuDifficultyHitObject lastObj = (OsuDifficultyHitObject)previousHandNotes[i + 1];

                double currHistoricalDecay = (history_time_max - (current.StartTime - currObj.StartTime)) / history_time_max; // scales note 0 to 1 from history to now

                currHistoricalDecay = Math.Min((double)(previousHandNotes.Count - i) / previousHandNotes.Count, currHistoricalDecay); // either we're limited by time or limited by object count.

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
                        if (previousHandNotes[i - 1].BaseObject is Slider) // bpm change is into slider, this is easy acc window
                            effectiveRatio *= 0.125;

                        if (previousHandNotes[i].BaseObject is Slider) // bpm change was from a slider, this is easier typically than circle -> circle
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
        private double calcAim(DifficultyHitObject current, DifficultyHitObject previous, DifficultyHitObject previousPrevious)
        {
            const double wide_angle_multiplier = 1.5;
            const double acute_angle_multiplier = 2.0;
            const double slider_multiplier = 1.5;
            const double velocity_change_multiplier = 0.75;

            if (current.BaseObject is Spinner || previousPrevious == null || previous.BaseObject is Spinner)
                return 0;

            var osuCurrObj = (OsuDifficultyHitObject)current;
            var osuLastObj = (OsuDifficultyHitObject)previous;
            var osuLastLastObj = (OsuDifficultyHitObject)previousPrevious;

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
        private double calcSpeed(DifficultyHitObject current, DifficultyHitObject previous)
        {
            const double single_spacing_threshold = 125;

            const double min_speed_bonus = 75; // ~200BPM
            const double speed_balancing_factor = 40;

            if (current.BaseObject is Spinner)
                return 0;

            // derive strainTime for calculation
            var osuCurrObj = (OsuDifficultyHitObject)current;
            var osuPrevObj = (OsuDifficultyHitObject)previous;

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
