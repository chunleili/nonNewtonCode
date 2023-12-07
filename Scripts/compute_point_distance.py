from scipy import spatial
import numpy as np
import meshio
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, default="spot_dense_scale.ply")
args = parser.parse_args()
filename = "data/models/" + args.filename

mesh = meshio.read(filename)
pos = mesh.points
tree = spatial.KDTree(mesh.points)

nearest_dist = np.zeros(shape=pos.shape[0])
nearest_indx = np.zeros(shape=pos.shape[0], dtype=np.int32)
for i in range(pos.shape[0]):
    nearest = tree.query(pos[i], k=2)
    nearest_dist[i] = nearest[0][1]
    nearest_indx[i] = nearest[1][1]
avg_distance = np.mean(nearest_dist)
max_distance = np.max(nearest_dist)
min_distance = np.min(nearest_dist)
print("total points: ", pos.shape[0])
print("max_distance: ", max_distance)
print("min_distance: ", min_distance)
print("avg_distance: ", avg_distance)
