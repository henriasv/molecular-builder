# Første linje
examples: brucite_in_water alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.png alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.png

.PHONY: brucite_in_water
brucite_in_water: brucite_in_water/brucite_in_water_orthogonal.png brucite_in_water/brucite_in_water_perspective.png

brucite_in_water/brucite_in_water_orthogonal.png brucite_in_water/brucite_in_water_perspective.png: brucite_in_water/brucite_in_water.py
	cd brucite_in_water && python brucite_in_water.py && cd ..

alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.png: alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.py
	cd alpha_quartz_cylinder_hole && python alpha_quartz_cylinder_hole.py && cd ..

alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.png: alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.py
	cd alpha_quartz_cylinder_hole_water && python alpha_quartz_cylinder_hole_water.py && cd ..

wedges_and_indenters/indenter_and_glass.png: wedges_and_indenters/indenter_over_silica.py
	cd wedges_and_indenters && python indenter_over_silica.py && cd ..

all: examples

