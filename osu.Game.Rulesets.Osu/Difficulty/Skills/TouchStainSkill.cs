// Copyright (c) ppy Pty Ltd <contact@ppy.sh>. Licensed under the MIT Licence.
// See the LICENCE file in the repository root for full licence text.

using System;
using osu.Game.Rulesets.Difficulty.Preprocessing;
using osu.Game.Rulesets.Difficulty.Utils;
using osu.Game.Rulesets.Mods;
using osu.Game.Rulesets.Objects;
using osu.Game.Rulesets.Osu.Difficulty.Preprocessing;

namespace osu.Game.Rulesets.Osu.Difficulty.Skills
{
    /// <summary>
    /// Represents the skill required to correctly aim and tap to the notes in a beatmap, using two hands on a touchscreen.
    /// </summary>
    public class TouchStrainSkill : OsuStrainSkill
    {
        protected override int HistoryLength => 60;
        private const int history_required = 32;
        private const int sequence_length = 4;

        // The actions utilized are bitmasked into an integer
        private const int l = 0;
        private const int r = 1;

        private const double overall_multiplier = 0.85; // Multiplier for final value

        // Constants for hand-coordination bonuses
        private const double coordination_aim_max_bonus = 0.65;
        private const double coordination_speed_max_bonus = 0.18; // Exact constant such that fully alternating between two hands rewards as much as streams
        private const double coordination_aim_decay = 5.6;
        private const double coordination_obstruction_decay = 1.8;
        private const double coordination_speed_decay = 9.5;
        private const float tap_x_obstruction_cap = 1.5f;

        // We keep an array of actions for convenience in code
        private readonly int[] actions =
        {
            l, r
        };

        private readonly double clockRate;
        private readonly bool calculatingAim;

        // Having the weighted strain of the most recent note stored globally to determine the initial strain of each strain section
        private double weightedAimStrain = 1;
        private double weightedSpeedStrain = 1;
        private double weightedRhythm = 1;

        // The aim, speed, and likelihood of each sequence are stored in an array with each index representing the bitmask of the sequence
        private readonly double[] sequenceAimStrain = new double[1 << sequence_length];
        private readonly double[] sequenceSpeedStrain = new double[1 << sequence_length];
        private readonly double[] likelihood = new double[1 << sequence_length]; // Should add up to 1.0
        private readonly long[] fullSequence = new long[1 << sequence_length]; // Contains full sequence up to HistoryLength notes to maintain approximate precision
        private readonly Aim aimSkill;
        private readonly Aim aimSkillNoSliders; // Obstruction bonuses should not stack with sliders
        private readonly Speed speedSkill;

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
            aimSkill = new Aim(mods, withSliders);
            aimSkillNoSliders = new Aim(mods, false);
            speedSkill = new Speed(mods, hitWindowGreat, true);
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
                return overall_multiplier * weightedAimStrain * decay(aimSkill.StrainDecayBase, timeDifference);

            return overall_multiplier * weightedRhythm * weightedSpeedStrain * decay(speedSkill.StrainDecayBase, timeDifference);
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
                double prevAim = sequenceAimStrain[i] * decay(aimSkill.StrainDecayBase, current.DeltaTime);
                double prevSpeed = sequenceSpeedStrain[i] * decay(speedSkill.StrainDecayBase, current.DeltaTime);
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
                    HitObject[] lastSame = new HitObject[history_required];
                    HitObject lastSwap = null;

                    // Rightmost bit that is NOT equal to "hand" corresponds to the index of the last swap
                    int lastSwapIndex = rightmostBitLocation(bitFixedExtension ^ invert);
                    if (lastSwapIndex < currentHistoryLength)
                        lastSwap = trueHistory[lastSwapIndex];

                    int mostRecentSame = 0;

                    while (mostRecentSame < history_required)
                    {
                        long leastSignificantBit = bitFixedExtension & -bitFixedExtension;
                        int lastSameIndex = rightmostBitLocation(bitFixedExtension);
                        // Only occurs if the current sequence of notes does not contain the current hand
                        if (lastSameIndex >= currentHistoryLength) break;

                        lastSame[mostRecentSame] = trueHistory[lastSameIndex];
                        mostRecentSame++;
                        // Subtracting LSB simulates lets rightmostBitLocation() identify the 2nd most recent note hit with the same hand, then the 3rd, etc...
                        bitFixedExtension -= leastSignificantBit;
                    }

                    if (mostRecentSame > 0)
                    {
                        // Calculating the individual strains
                        var simulatedCurrent = new OsuDifficultyHitObject(current.BaseObject, lastSame[1], lastSame[0], clockRate);
                        var simulatedPrevious = new ReverseQueue<DifficultyHitObject>(history_required);
                        const int min_required = 3; // Minimum DifficultyHitObjects required for aim/speed calculation, the rest will be "lazily" generated, only consisting of deltaTimes

                        for (int j = mostRecentSame - 2; j > min_required; j--)
                            simulatedPrevious.Enqueue(new OsuDifficultyHitObject(lastSame[j], j + 2 == mostRecentSame ? null : lastSame[j + 2], lastSame[j + 1], clockRate, false));

                        for (int j = Math.Min(min_required, mostRecentSame - 2); j >= 0; j--)
                            simulatedPrevious.Enqueue(new OsuDifficultyHitObject(lastSame[j], j + 2 == mostRecentSame ? null : lastSame[j + 2], lastSame[j + 1], clockRate));

                        // Adding a hand-coordination bonus based on time since the last hand "swap"
                        double aimBonus = 1.0;
                        double speedBonus = 1.0;
                        double obstructionBonus = 1.0;

                        if (lastSwap != null)
                        {
                            var simulatedSwap = new OsuDifficultyHitObject(current.BaseObject, lastSame[0], lastSwap, clockRate);
                            double ratio = simulatedCurrent.DeltaTime / (simulatedCurrent.DeltaTime + simulatedSwap.DeltaTime);

                            // Add to obstruction bonus if hands cross over each other horizontally, assuming the right hand is most comfortable being to the right of the left hand
                            double handAdjust = hand == l ? 1 : -1;
                            double horizontalOverlap = simulatedSwap.HorizontalDisplacementFactor * handAdjust;
                            horizontalOverlap = 1.0 / (1.0 + Math.Pow(Math.E, 2.5 - 6.0 * horizontalOverlap));

                            double rawObstruction = 2.0 * simulatedSwap.ObstructionFactor;
                            obstructionBonus += (horizontalOverlap + rawObstruction) * coordinationFactor(ratio, coordination_obstruction_decay);

                            aimBonus += coordination_aim_max_bonus * coordinationFactor(ratio, coordination_aim_decay);
                            speedBonus += coordination_speed_max_bonus * coordinationFactor(ratio, coordination_speed_decay);
                        }

                        // Do not apply obstruction bonuses to slider travel time.
                        double rawAimSkill = aimSkill.StrainValueOf(simulatedCurrent, simulatedPrevious);
                        double rawMovement = aimSkillNoSliders.StrainValueOf(simulatedCurrent, simulatedPrevious);
                        simulatedCurrent.UpdateDistances((float)obstructionBonus);

                        for (int j = Math.Min(min_required, mostRecentSame - 2); j >= 0; j--)
                        {
                            var previous = (OsuDifficultyHitObject)simulatedPrevious[j];
                            previous.UpdateDistances((float)obstructionBonus);
                        }

                        rawMovement = aimSkillNoSliders.StrainValueOf(simulatedCurrent, simulatedPrevious) - rawMovement;

                        // Set a hard limit on obstruction bonuses to take tapX into consideration.
                        if (obstructionBonus >= tap_x_obstruction_cap)
                        {
                            simulatedCurrent.UpdateDistances(tap_x_obstruction_cap);

                            for (int j = Math.Min(min_required, mostRecentSame - 2); j >= 0; j--)
                            {
                                var previous = (OsuDifficultyHitObject)simulatedPrevious[j];
                                previous.UpdateDistances(tap_x_obstruction_cap);
                            }
                        }

                        double rawSpeedSkill = speedSkill.StrainValueOf(simulatedCurrent, simulatedPrevious);

                        rawAim[hand] = aimSkill.SkillMultiplier * aimBonus * (rawMovement + rawAimSkill);
                        rawSpeed[hand] = speedSkill.SkillMultiplier * speedBonus * rawSpeedSkill;
                        rawRhythm[hand] = speedSkill.CalculateRhythmBonus(simulatedCurrent, simulatedPrevious);
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

            weightedRhythm = Math.Sqrt(weightedRhythm * speedSkill.CalculateRhythmBonus(current, Previous));
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

        private double coordinationFactor(double swapRatio, double decay) => Math.Min(1, Math.Pow(2.0, decay) * Math.Pow(swapRatio, decay));

        private double curveProbability(double chance) => chance <= 0.5 ? Math.Pow(2 * chance, 4) / 2 : 1 - Math.Pow(2 * (1 - chance), 4) / 2;

        private double decay(double decayBase, double ms) => Math.Pow(decayBase, ms / 1000);

        private double curveStrain(double aimStrain, double speedStrain) => Math.Pow(Math.Pow(aimStrain, 3.0 / 2.0) + Math.Pow(speedStrain, 3.0 / 2.0), 2.0 / 3.0);

        private int appendAction(int oldSequence, int toAppend) => ((oldSequence << 1) | toAppend) & ((1 << sequence_length) - 1);

        private long appendHistory(long oldSequence, long toAppend) => ((oldSequence << 1) | toAppend) & ((1 << HistoryLength) - 1);

        private int rightmostBitLocation(long value) => fastLog2((ulong)(value & -value));

        private int fastLog2(ulong value) => logTable[(int)((value * 0x07EDD5E59A4E28C2) >> 58)];
    }
}
