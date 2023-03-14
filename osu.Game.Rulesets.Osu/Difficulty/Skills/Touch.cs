// Copyright (c) ppy Pty Ltd <contact@ppy.sh>. Licensed under the MIT Licence.
// See the LICENCE file in the repository root for full licence text.

using System;
using System.Collections.Generic;
using System.Linq;
using osu.Game.Rulesets.Difficulty.Preprocessing;
using osu.Game.Rulesets.Mods;
using osu.Game.Rulesets.Objects;
using osu.Game.Rulesets.Osu.Difficulty.Evaluators;
using osu.Game.Rulesets.Osu.Difficulty.Preprocessing;

namespace osu.Game.Rulesets.Osu.Difficulty.Skills
{
    /// <summary>
    /// Represents the skill required to correctly aim and tap to the notes in a beatmap, using two hands on a touchscreen.
    /// </summary>
    public class Touch : OsuStrainSkill
    {
        private const int l = 0;
        private const int r = 1;
        private const int max_possibilities = 15;

        // Constants from aim and speed for standard gameplay
        private const double aim_skill_multiplier = 23.55;
        private const double aim_strain_decay_base = 0.15;
        private const double speed_skill_multiplier = 1375;
        private const double speed_strain_decay_base = 0.3;
        private class Possibility
        {
            private const int difficulty_objects_required = 3;
            private const int beatmap_objects_required = 2;
            public double Probability;
            public double Aim { get; private set; }
            public double Speed { get; private set; }
            private readonly double clockRate;
            private readonly bool withSliders;
            private readonly List<DifficultyHitObject>[] objects = new List<DifficultyHitObject>[2];
            private readonly List<HitObject>[] beatmap = new List<HitObject>[2];
            private int lastHand;
            private int otherHand(int currentHand) => 1 - currentHand;
            public Possibility(double clockRate, bool withSliders, DifficultyHitObject firstObject)
            {
                this.clockRate = clockRate;
                this.withSliders = withSliders;
                this.Probability = 1;
                objects[l] = new List<DifficultyHitObject>();
                objects[r] = new List<DifficultyHitObject>();

                beatmap[l] = new List<HitObject>();
                beatmap[r] = new List<HitObject>();

                // Automatically assume the first note of a beatmap is hit with the left hand and the second note is hit with the right.
                beatmap[l].Add(firstObject.LastObject);
                beatmap[r].Add(firstObject.BaseObject);
                lastHand = r;
            }
            public Possibility(Possibility copy)
            {
                this.clockRate = copy.clockRate;
                this.withSliders = copy.withSliders;
                this.Probability = copy.Probability;
                this.lastHand = copy.lastHand;
                this.Aim = copy.Aim;
                this.Speed = copy.Speed;

                objects[l] = new List<DifficultyHitObject>(copy.objects[l]);
                objects[r] = new List<DifficultyHitObject>(copy.objects[r]);

                beatmap[l] = new List<HitObject>(copy.beatmap[l]);
                beatmap[r] = new List<HitObject>(copy.beatmap[r]);
            }
            public void UpdateStrainValue(DifficultyHitObject current, int currentHand)
            {
                double strainTime = ((OsuDifficultyHitObject)current).StrainTime;
                Aim *= strainDecay(aim_strain_decay_base, strainTime);
                Speed *= strainDecay(aim_strain_decay_base, strainTime);

                double aimIfCurrentHand = aimStrainValueIf(current, currentHand);
                double speedIfCurrentHand = speedStrainValueIf(current, currentHand);

                Aim += aimIfCurrentHand;
                Speed += speedIfCurrentHand;

                double aimIfOtherHand = aimStrainValueIf(current, otherHand(currentHand));
                double speedIfOtherHand = speedStrainValueIf(current, otherHand(currentHand));

                double strainIfCurrentHand = totalStrain(aimIfCurrentHand, speedIfCurrentHand);
                double strainIfOtherHand = totalStrain(aimIfOtherHand, speedIfOtherHand);

                double probabilityCurrentHand = strainIfCurrentHand > 0 ? strainIfOtherHand / (strainIfCurrentHand + strainIfOtherHand) : 1;

                // Curve probability to be closer to either 0 or 1.
                if (probabilityCurrentHand <= 0.5)
                    probabilityCurrentHand = 2 * Math.Pow(probabilityCurrentHand, 2);
                else
                    probabilityCurrentHand = 1 - 2 * Math.Pow(1 - probabilityCurrentHand, 2);

                Probability *= probabilityCurrentHand;
            }
            public void UpdateHistory(DifficultyHitObject current, int currentHand)
            {
                beatmap[currentHand].Add(current.BaseObject);
                objects[currentHand].Add(getSimulatedObject(current, currentHand));

                if (beatmap[currentHand].Count > beatmap_objects_required)
                    beatmap[currentHand].RemoveRange(0, beatmap[currentHand].Count - beatmap_objects_required);

                if (objects[currentHand].Count > difficulty_objects_required)
                    objects[currentHand].RemoveRange(0, objects[currentHand].Count - difficulty_objects_required);

                lastHand = currentHand;
            }
            private OsuDifficultyHitObject getSimulatedObject(DifficultyHitObject current, int currentHand)
            {
                // A simulated difficulty object is created for angle calculation.
                int index = beatmap[currentHand].Count;
                var lastLast = index > 1 ? beatmap[currentHand][index - 2] : null;
                return new OsuDifficultyHitObject(current.BaseObject, beatmap[currentHand][index - 1], lastLast, clockRate, objects[currentHand], objects[currentHand].Count);
            }
            private OsuDifficultyHitObject getSimulatedSwapObject(DifficultyHitObject current, int currentHand)
            {
                int i1 = beatmap[currentHand].Count;
                int i2 = beatmap[otherHand(currentHand)].Count;
                var last = beatmap[otherHand(currentHand)][i2 - 1];
                var lastLast = beatmap[currentHand][i1 - 1];
                return new OsuDifficultyHitObject(current.BaseObject, last, lastLast, clockRate, objects[currentHand], objects[currentHand].Count);
            }
            private double aimStrainValueIf(DifficultyHitObject current, int currentHand)
            {
                double obstructionBonus = 1;
                // Add a bonus for the hand co-ordination required to aim with both hands.
                if (currentHand != lastHand)
                {
                    obstructionBonus += 1.25;

                    // Add an obstrution bonus if the most recent instance of the "other hand" is in between the current object and the previous object with the actual hand
                    var simulatedSwap = getSimulatedSwapObject(current, currentHand);
                    double? angle = simulatedSwap.Angle;
                    if (angle != null)
                        obstructionBonus += 1.25 / (1 + Math.Pow(Math.E, -(angle.Value * 180 / Math.PI - 108) / 9));
                }

                return AimEvaluator.EvaluateDifficultyOf(getSimulatedObject(current, currentHand), withSliders) * obstructionBonus * aim_skill_multiplier;
            }
            private double speedStrainValueIf(DifficultyHitObject current, int currentHand)
            {
                return SpeedEvaluator.EvaluateDifficultyOf(current, false) * speed_skill_multiplier;
            }
            private double totalStrain(double aimStrain, double speedStrain) => Math.Pow(Math.Pow(aimStrain, 3.0 / 2.0) + Math.Pow(speedStrain, 3.0 / 2.0), 2.0 / 3.0);
        }

        private readonly bool calculatingAim;
        private readonly double clockRate;
        private readonly bool withSliders;

        private double weightedAim;
        private double weightedSpeed;
        private double rhythm;

        private readonly List<double> objectStrains = new List<double>();
        private List<Possibility> possibilities = new List<Possibility>();

        public Touch(Mod[] mods, double clockRate, bool withSliders, bool calculatingAim)
            : base(mods)
        {
            this.clockRate = clockRate;
            this.withSliders = withSliders;
            this.calculatingAim = calculatingAim;
        }

        protected override double CalculateInitialStrain(double time, DifficultyHitObject current)
        {
            double timeDifference = time - current.Previous(0).StartTime;
            if (calculatingAim)
                return weightedAim * strainDecay(aim_strain_decay_base, timeDifference);

            return rhythm * weightedSpeed * strainDecay(speed_strain_decay_base, timeDifference);
        }

        protected override double StrainValueAt(DifficultyHitObject current)
        {
            weightedAim = 0;
            weightedSpeed = 0;
            rhythm = RhythmEvaluator.EvaluateDifficultyOf(current);

            if (possibilities.Count == 0)
            {
                var initialPossibility = new Possibility(clockRate, withSliders, current);
                possibilities.Add(initialPossibility);
                return 0;
            }

            List<Possibility> newPossibilities = new List<Possibility>();

            foreach (var cur in possibilities)
            {
                var left = new Possibility(cur);
                var right = new Possibility(cur);

                left.UpdateStrainValue(current, l);
                right.UpdateStrainValue(current, r);

                left.UpdateHistory(current, l);
                right.UpdateHistory(current, r);

                newPossibilities.Add(left);
                newPossibilities.Add(right);
            }

            // Only keep the most probable possibilities.
            newPossibilities = newPossibilities.OrderByDescending(x => x.Probability).ToList();
            possibilities = new List<Possibility>();

            int entries = Math.Min(newPossibilities.Count, max_possibilities);
            double totalMostProbable = 0;

            for (int i = 0; i < entries; i++)
            {
                totalMostProbable += newPossibilities[i].Probability;
                possibilities.Add(newPossibilities[i]);
            }

            // Make sure total probability sums up to 1.
            foreach (var cur in possibilities)
                cur.Probability = totalMostProbable > 0 ? cur.Probability / totalMostProbable : 1 / entries;

            foreach (var x in possibilities)
            {
                weightedAim += x.Aim * x.Probability;
                weightedSpeed += x.Speed * x.Probability;
            }

            if (calculatingAim)
                return weightedAim;

            double totalStrain = rhythm * weightedSpeed;

            objectStrains.Add(totalStrain);

            return totalStrain;
        }
        public double RelevantNoteCount()
        {
            if (objectStrains.Count == 0)
                return 0;

            double maxStrain = objectStrains.Max();

            if (maxStrain == 0)
                return 0;

            return objectStrains.Sum(strain => 1.0 / (1.0 + Math.Exp(-(strain / maxStrain * 12.0 - 6.0))));
        }
        private static double strainDecay(double decayBase, double ms) => Math.Pow(decayBase, ms / 1000);
    }
}