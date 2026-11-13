import pickle
import itertools

with open("shortest_paths_cache.pkl", "rb") as f:
    cache = pickle.load(f)

print(type(cache))
print("Total nodes:", len(cache))

for node, distances in itertools.islice(cache.items(), 10):
    print("\nNode:", node)
    print("First 20 shortest-path lengths:")

    for target, path_length in itertools.islice(distances.items(), 20):
        print(f"  {target}: {path_length}")