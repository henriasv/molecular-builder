# Første linje
examples: cracks procedural_noise sodium_chlorate_hydrocarbon_blob_water brucite_in_water alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.png alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.png wedges_and_indenters/indenter_over_silica.png silicon_carbide

brucite_in_water: brucite_in_water/brucite_in_water_orthogonal.png brucite_in_water/brucite_in_water_perspective.png

brucite_in_water/brucite_in_water_orthogonal.png brucite_in_water/brucite_in_water_perspective.png: brucite_in_water/brucite_in_water.py
	cd brucite_in_water && python brucite_in_water.py && cd ..

alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.png: alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.py
	cd alpha_quartz_cylinder_hole && python alpha_quartz_cylinder_hole.py && cd ..

alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.png: alpha_quartz_cylinder_hole_water/alpha_quartz_cylinder_hole_water.py
	cd alpha_quartz_cylinder_hole_water && python alpha_quartz_cylinder_hole_water.py && cd ..

wedges_and_indenters/indenter_over_silica.png: wedges_and_indenters/indenter_over_silica.py
	cd wedges_and_indenters && python indenter_over_silica.py && cd ..

sodium_chlorate_hydrocarbon_blob_water: sodium_chlorate_hydrocarbon_blob_water/sodium_chlorate_hydrocarbon_blob_water.png

sodium_chlorate_hydrocarbon_blob_water/sodium_chlorate_hydrocarbon_blob_water.png: sodium_chlorate_hydrocarbon_blob_water/sodium_chlorate_hydrocarbon_blob_water.py
	cd sodium_chlorate_hydrocarbon_blob_water && python sodium_chlorate_hydrocarbon_blob_water.py && cd ..

silicon_carbide: silicon_carbide_nanoparticles/box.png silicon_carbide_nanoparticles/cylinder.png silicon_carbide_nanoparticles/sphere.png silicon_carbide_nanoparticles/dodecahedron.png silicon_carbide_nanoparticles/octahedron.png

silicon_carbide_nanoparticles/box.png: silicon_carbide_nanoparticles/box.py
	cd silicon_carbide_nanoparticles && python box.py && cd ..

silicon_carbide_nanoparticles/sphere.png: silicon_carbide_nanoparticles/sphere.py
	cd silicon_carbide_nanoparticles && python sphere.py && cd ..

silicon_carbide_nanoparticles/cylinder.png: silicon_carbide_nanoparticles/cylinder.py
	cd silicon_carbide_nanoparticles && python cylinder.py && cd ..

silicon_carbide_nanoparticles/octahedron.png: silicon_carbide_nanoparticles/octahedron.py
	cd silicon_carbide_nanoparticles && python octahedron.py && cd ..

silicon_carbide_nanoparticles/dodecahedron.png: silicon_carbide_nanoparticles/dodecahedron.py
	cd silicon_carbide_nanoparticles && python dodecahedron.py && cd ..

procedural_noise/block.png: procedural_noise/generate.py
	cd procedural_noise && python generate.py && cd ..

procedural_noise/basic.png: procedural_noise/basic.py
	cd procedural_noise && python basic.py && cd ..

procedural_noise/parameters.png: procedural_noise/parameters.py
	cd procedural_noise && python parameters.py && cd ..

procedural_noise/two_level.png: procedural_noise/two_level.py
	cd procedural_noise && python two_level.py && cd ..

procedural_noise/periodic.png: procedural_noise/periodic.py
	cd procedural_noise && python periodic.py && cd ..

procedural_noise/periodic2.png: procedural_noise/periodic2.py
	cd procedural_noise && python periodic2.py && cd ..

procedural_noise/screwed.png: procedural_noise/screwed.py
	cd procedural_noise && python screwed.py && cd ..

procedural_noise/add_function.png: procedural_noise/add_function.py
	cd procedural_noise && python add_function.py && cd ..

cracks: cracks/alpha_quartz_with_crack.png

cracks/alpha_quartz_with_crack.png: cracks/crack_in_quartz.py
	cd cracks && python crack_in_quartz.py && cd ..

water_bulk/box.png: water_bulk/box.py
	cd water_bulk && python box.py && cd ..

all: examples
