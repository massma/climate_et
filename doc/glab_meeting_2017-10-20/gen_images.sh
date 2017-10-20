rm -rf ./*.png
convert ../paper/fig04.pdf -gravity South -crop 100x50%+0+0 ./scaling.png
convert ../paper/fig05.pdf -crop 100x50%+0+0 ./sign.png
convert ../paper/fig06.png -gravity NorthEast -crop 25x20%+0+0 ./temp_1.png
convert ../paper/fig06b.png -gravity NorthEast -crop 25x20%+0+0 ./temp_2.png
convert +append temp*.png cro_ideal.png
rm -rf ./temp*.png
convert ../paper/fig06.png -crop 25x20%+0+0 ./temp_1.png
convert ../paper/fig06b.png -crop 25x20%+0+0 ./temp_2.png
convert +append temp*.png cro.png
rm -rf ./temp*.png

# idealized plots
convert ../paper/fig06.png -gravity East -crop 25x20%+0+0 ./temp_1.png
convert ../paper/fig06b.png -gravity East -crop 25x20%+0+0 ./temp_2.png
convert +append temp*.png dbf_ideal.png
rm -rf ./temp*.png
convert ../paper/fig06.png -gravity West -crop 25x20%+0+0 ./temp_1.png
convert ../paper/fig06b.png -gravity West -crop 25x20%+0+0 ./temp_2.png
convert +append temp*.png dbf.png
rm -rf ./temp*.png


# idealized plots
convert ../paper/fig06.png -crop 640x480+0+480 ./temp_1.png
convert ../paper/fig06b.png -crop 640x480+0+480 ./temp_2.png
convert +append temp*.png csh.png
rm -rf ./temp*.png
convert ../paper/fig06.png -gravity NorthEast -crop 640x480+0+480 ./temp_1.png
convert ../paper/fig06b.png -gravity NorthEast -crop 640x480+0+480 ./temp_2.png
convert +append temp*.png csh_ideal.png
rm -rf ./temp*.png


convert ../paper/fig06.png -crop 640x480+0+1440 ./temp_1.png
convert ../paper/fig06b.png -crop 640x480+0+1440 ./temp_2.png
convert +append temp*.png enf.png
rm -rf ./temp*.png
convert ../paper/fig06.png -gravity NorthEast -crop 640x480+0+1440 ./temp_1.png
convert ../paper/fig06b.png -gravity NorthEast -crop 640x480+0+1440 ./temp_2.png
convert +append temp*.png enf_ideal.png
rm -rf ./temp*.png

convert ../paper/fig06.png -crop 640x480+0+1920 ./temp_1.png
convert ../paper/fig06b.png -crop 640x480+0+1920 ./temp_2.png
convert +append temp*.png gra.png
rm -rf ./temp*.png
convert ../paper/fig06.png -gravity NorthEast -crop 640x480+0+1920 ./temp_1.png
convert ../paper/fig06b.png -gravity NorthEast -crop 640x480+0+1920 ./temp_2.png
convert +append temp*.png gra_ideal.png
rm -rf ./temp*.png
