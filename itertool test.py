from itertools import combinations

# Individual cuts (assuming these are pandas DataFrame filters)
individualCuts = [cutDD, cutOM, cutPLA, cutBAO, cutCMB, cutLZ, cutEW]
cutList = []

# Generate all non-empty subsets while maintaining order
for r in range(1, len(individualCuts) + 1):  # r is the size of the combination
    for combo in combinations(individualCuts, r):
        cutList.append(" & ".join(combo))

# Print the result
for cut in cutList:
    print(cut)