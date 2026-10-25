from dnachisel import *


protein_seq = "MNLEGGGRGGEFGMSAVSCGNGKLRQWLIDQIDSGKYPGLVWENEEKSIFRIPWKHAGKQDYNREEDAALFKAWALFKGKFREGIDKPDPPTWKTRLRCALNKSNDFEELVERSQLDISDPYKVYRIVPEGAKKGAKQLTLEDPQMSMSHPYTMTTPYPSLPAQQVHNYMMPPLDRSWRDYVPDQPHPEIPYQCPMTFGPRGHHWQGPACENGCQVTGTFYACAPPESQAPGVPTEPSIRSAEALAFSDCRLHICLYYREILVKELTTSSPEGCRISHGHTYDASNLDQVLFPYPEDNGQRKNIEKLLSHLERGVVLWMAPDGLYAKRLCQSRIYWDGPLALCNDRPNKLERDQTCKLFDTQQFLSELQAFAHHGRSLPRFQVTLCFGEEFPDPQRQRKLITAHVEPLLARQLYYFAQQNSGHFLRGYDLPEHISNPEDYHRSIRHSSIQE"

some_sequence = reverse_translate(protein_seq)

problem = DnaOptimizationProblem(
    sequence=some_sequence,
    constraints=[
        EnforceTranslation(),
        AvoidHairpins(),

        AvoidPattern("AAAAA"),
        AvoidPattern("TTTTT"),
        AvoidPattern("GGGGG"),
        AvoidPattern("CCCCC"),

        AvoidRareCodons(min_frequency=0.1,species='h_sapiens'),
        EnforceGCContent(mini=0.3, maxi=0.7, window=50),
    ],
    objectives=[CodonOptimize(species='h_sapiens')]
)

# SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

problem.resolve_constraints()
problem.optimize()

# PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

print(problem.constraints_text_summary())
print(problem.objectives_text_summary())
print(problem.sequence)
