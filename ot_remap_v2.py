#!/usr/bin/env python2
# ot_remap_v2.py -- creates a remapped model (quickly) from a RELION subtomogram average
# Created 20210329 (Jon Chen)
# Revised 20210503 (LG) Added user-friendly inputs.
# Revised 20220208 (LG) Add fixes for python3 and EMAN2.91 (uses Python 3.7)
# Revised 20220822 (JC) Minor fix to rotation handling of translational offsets (affected line commented "fix001")

import mrcfile, numpy, copy, math
from EMAN2 import *

#---------- BEGIN User inputs ------------------------------
# stackoverflow.com/questions/2194163/python-empty-argument
try:
   avg_file  = sys.argv[1]
   tomo_size = sys.argv[2]
   star_file = sys.argv[3]
   out_file  = sys.argv[4]
except IndexError:
   print("================================================================")
   print("Usage:>>   ot_remap_v2.py [Average] [Tomogram_size] [Starfile] [Remapped]")
   print("Example:   ot_remap_v2.py run_class001.mrc 2048,2048,500 run_data.star remap_001.mrc")
   print("----------------------------------------------------------------")
   print("Make sure EMAN2's binaries are in the PATH using a command like:")
   print("export PATH=/mnt/prog/EMAN2.91/bin:$PATH")
   print("----------------------------------------------------------------")
   print("Average:   The class average you want to remap")
   print("Tomogram_size:  Output tomogram size (XYZ, comma-separated)")
   print("Starfile:  _data.star file that contains the refined coordinates & orientations")
   print("Remapped:  Name of the remapped .mrc file")
   print("----------------------------------------------------------------")
   print("Output:>> [Remapped]")
   print("To invert the remapped tomogram, run 'bimg -invert [Remapped] negative.mrc'")
   sys.exit()
#----------- END User inputs -------------------------------


def get_star_index(filepath):
	col_index = {}
	with open(filepath, "r") as f:
		for i in range(10):
			line = f.readline().strip()
			if line == "loop_":
				break
			elif i == 9:
				print("Could not detect loop_ line.. exiting")
				exit()
		while 1:
			line = f.readline().strip()
			if line[0] != "_":
				break
			else:
				x = line[1:].split()
				col_index[x[0]] = int(x[1].replace("#", "")) - 1
	return col_index

def get_star_data(filepath):
	data_list = []
	with open(filepath, "r") as f:
		for line in f:
			if line.strip() == "loop_":
				break
		for line in f:
			line_strip = line.strip()
			try:
				line_strip[0]
			except IndexError:
				continue
			if line_strip[0] == "_":
				continue
			else:
				row = " ".join(line_strip.split()).split(" ")
				data_list.append(row)
	return data_list

def get_particle_data(star_file):
	col_index = get_star_index(star_file)
	data_list = get_star_data(star_file)
	x_index = col_index["rlnCoordinateX"]
	y_index = col_index["rlnCoordinateY"]
	z_index = col_index["rlnCoordinateZ"]
	x_offset_index = col_index["rlnOriginX"]
	y_offset_index = col_index["rlnOriginY"]
	z_offset_index = col_index["rlnOriginZ"]
	rot_index = col_index["rlnAngleRot"]
	tilt_index = col_index["rlnAngleTilt"]
	psi_index = col_index["rlnAnglePsi"]
	particle_data = []
	for row in data_list:
		refined_x = float(row[x_index]) - float(row[x_offset_index])
		refined_y = float(row[y_index]) - float(row[y_offset_index])
		refined_z = float(row[z_index]) - float(row[z_offset_index])
		rot_angle = float(row[rot_index])
		tilt_angle = float(row[tilt_index])
		psi_angle = float(row[psi_index])
		particle_data.append([refined_x, refined_y, refined_z, rot_angle, tilt_angle, psi_angle])
	return particle_data

def rotation_matrix(axis, theta):
    """
    # modified version from https://stackoverflow.com/a/6802723
    Return the rotation matrix associated with CLOCKWISE rotation about
    the given axis by theta radians.
    """
    axis = numpy.asarray(axis)
    axis = axis / math.sqrt(numpy.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return numpy.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def zyz_rot(point, z1, y2, z3):
    # z3 around z-axis
    if z3 != 0:
        theta = numpy.radians(z3)
        axis = [0, 0, 1]
        point_rot = numpy.dot(rotation_matrix(axis, theta), point)
        point_rot = [round(num, 6) for num in point_rot]
        point = point_rot
    #y2 around y-axis
    if y2 != 0:
        theta = numpy.radians(y2)
        axis = [0, 1, 0]
        point_rot = numpy.dot(rotation_matrix(axis, theta), point)
        point_rot = [round(num, 6) for num in point_rot]
        point = point_rot
    #z1 around z-axis
    if z1 != 0:
        theta = numpy.radians(z1)
        axis = [0, 0, 1]
        point_rot = numpy.dot(rotation_matrix(axis, theta), point)
        point_rot = [round(num, 6) for num in point_rot]
        point = point_rot
    return point

def nearest_half(number):
	number_round = round(number*2)/2
	if number_round%1 != 0.5:
		if number_round > number:
			number_round -= 0.5
		elif number_round < number:
			number_round += 0.5
	return number_round

def main(avg_file, tomo_file, star_file, out_file):
	# Parse tomo_size argument
	tomo_xyz_size = tomo_size.split(",")
	tomo_xyz_size = [int(val) for val in tomo_xyz_size]

	# Get dimensions of subtomo average file
	avg_xyz_size = [0, 0, 0]
	with mrcfile.mmap(avg_file, mode='r') as tomo:
		avg_xyz_size[0] = tomo.header.nx.item()
		avg_xyz_size[1] = tomo.header.ny.item()
		avg_xyz_size[2] = tomo.header.nz.item()

	# Work out the offset coordinate of subtomo average center 
	# with respect to EMAN2's rotation center:
	boxsize = avg_xyz_size[0]
	avg_center = [0, 0, 0]
	for i in range(len(avg_xyz_size)):
		if avg_xyz_size[i]%2 == 0:
			avg_center[i] = -0.5
		elif avg_xyz_size[i]%2 == 1:
			avg_center[i] = 0

	# Get particle coordinates and orientation from star file
	particle_data = get_particle_data(star_file)

	# Read subtomo average file using EMAN2 functions
	avg = EMData()
	avg.read_image(avg_file)

	# Create output mrc file as memory-mapped object using mrcfile
	tomo = mrcfile.new_mmap(out_file, shape=(tomo_xyz_size[2], tomo_xyz_size[1], tomo_xyz_size[0]), mrc_mode=2)

	# Start remapping
	k = 0
	for row in particle_data:
		k += 1
		xyz_offset = [0, 0, 0]
		avg_xform = copy.copy(avg)
		coord = [row[0], row[1], row[2]]
		rot = [row[3], row[4], row[5]]

		# Apply rotation to subtomo average array
		t = Transform({"type":"spider", "phi":rot[0], "theta":rot[1], "psi":rot[2]})
		avg_xform.transform(t)

		# Work out xyz offset for rotated subtomo average array
		## If EMAN2 rotation center != true volume center, include xyz offset
		avg_center_xform = zyz_rot(avg_center, rot[2], rot[1], rot[0])  #fix001 - changed to inverse rotation (20220822)
		xyz_offset = [(avg_center[i] - avg_center_xform[i]) for i in range(len(avg_center))]
		## Round RELION coordinates to the closest 0.5, add difference to xyz offset
		coord_round = [nearest_half(coord[i]) for i in range(len(coord))]
		for i in range(len(xyz_offset)):
			xyz_offset[i] = xyz_offset[i] + coord[i] - coord_round[i]
		## Apply xyz offset to subtomo average array
		t = Transform({"tx":xyz_offset[0], "ty":xyz_offset[1], "tz":xyz_offset[2]})
		avg_xform.transform(t)

		# Convert subtomo average EMAN2 object to numpy array
		avg_arr = EMNumPy.em2numpy(avg_xform)
		avg_arr = avg_arr.astype("float32")
		# Work out coordinates in output tomo to write particle array
		# Will ignore writing parts of subtomo array that are outside of tomo volume
		# Currently only supports subtomo average volumes with same xyz-dimensions and even size
		corner_coord_min = [int(coord_round[i]-0.5-boxsize/2+1) for i in range(len(coord_round))]
		corner_coord_max = [corner_coord_min[i]+boxsize for i in range(len(corner_coord_min))]
		crop_min = [0, 0, 0]
		crop_max = [0, 0, 0]
		for i in range(len(corner_coord_min)):
			if corner_coord_min[i] < 0:
				crop_min[i] = 0 - corner_coord_min[i]
		for i in range(len(corner_coord_max)):
			if corner_coord_max[i] > tomo_xyz_size[i]:
				crop_max[i] = corner_coord_max[i] - tomo_xyz_size[i]
		a = [corner_coord_min[i] + crop_min[i] for i in range(len(corner_coord_min))]
		b = [corner_coord_max[i] - crop_max[i] for i in range(len(corner_coord_max))]
		c = [0 + crop_min[i] for i in range(len(corner_coord_min))]
		d = [boxsize - crop_max[i] for i in range(len(corner_coord_max))]
		
		# Write final subtomo volume array to out tomo
		## If multiple volumes have overlapping voxels in output tomo,
		## the largest value will be retained:
		mask = numpy.ma.masked_not_equal(tomo.data[a[2]:b[2], a[1]:b[1], a[0]:b[0]], 0, copy=True)
		masked_avg_arr = numpy.ma.masked_where(numpy.ma.getmask(mask), avg_arr[c[2]:d[2], c[1]:d[1], c[0]:d[0]])
		masked_avg_arr = numpy.ma.filled(masked_avg_arr, fill_value=9999)
		temp_arr = numpy.minimum.reduce([masked_avg_arr, tomo.data[a[2]:b[2], a[1]:b[1], a[0]:b[0]]])
		tomo.data[a[2]:b[2], a[1]:b[1], a[0]:b[0]] = numpy.maximum.reduce([avg_arr[c[2]:d[2], c[1]:d[1], c[0]:d[0]], temp_arr])

	print("Remapped %s particles!"%k)
	print("Updating tomo stats..")
	tomo.update_header_stats()
	print("Closing tomo..")
	tomo.close()

if __name__ == "__main__":
	avg_mrc = sys.argv[1]
	tomo_mrc = sys.argv[2]
	avg_star = sys.argv[3]
	out_mrc = sys.argv[4]
	main(avg_mrc, tomo_mrc, avg_star, out_mrc)
