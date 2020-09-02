# Første linje
examples: sodium_chlorate_hydrocarbon_blob_water brucite_in_water alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.png alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.png wedges_and_indenters/indenter_and_glass.png silicon_carbide

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

sodium_chlorate_hydrocarbon_blob_water: sodium_chlorate_hydrocarbon_blob_water/sodium_chlorate_hydrocarbon_blob_water.png

sodium_chlorate_hydrocarbon_blob_water/sodium_chlorate_hydrocarbon_blob_water.png: sodium_chlorate_hydrocarbon_blob_water/sodium_chlorate_hydrocarbon_blob_water.py
	cd sodium_chlorate_hydrocarbon_blob_water && python sodium_chlorate_hydrocarbon_blob_water.py && cd ..

silicon_carbide: silicon_carbide_nanoparticles/box.png silicon_carbide_nanoparticles/cylinder.png silicon_carbide_nanoparticles/sphere.png

silicon_carbide_nanoparticles/box.png: 
	cd silicon_carbide_nanoparticles && python box.py && cd ..

silicon_carbide_nanoparticles/sphere.png:
	cd silicon_carbide_nanoparticles && python sphere.py && cd ..

silicon_carbide_nanoparticles/cylinder.png:
	cd silicon_carbide_nanoparticles && python cylinder.py && cd ..

all: examples

