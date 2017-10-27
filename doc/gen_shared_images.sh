rm -rf ./shared_figs/*.png
convert ./paper/fig04.pdf -gravity South -crop 100x50%+0+0 ./shared_figs/scaling.png
convert ./paper/fig05.pdf -crop 100x48%+0+0 ./shared_figs/sign.png
pdfjam --keepinfo --trim "0mm 160mm 0mm 0mm" --clip true --suffix "cropped" ./paper/fig05.pdf
mv ./fig05-cropped.pdf ./shared_figs/sign.pdf
pdfcrop --margins '5 5 5 10' ./shared_figs/sign.pdf ./shared_figs/sign.pdf
convert ./paper/fig06.png -gravity NorthEast -crop 25x20%+0+0 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -gravity NorthEast -crop 25x20%+0+0 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/cro_ideal.png
rm -rf ./shared_figs/temp*.png
convert ./paper/fig06.png -crop 25x20%+0+0 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -crop 25x20%+0+0 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/cro.png
rm -rf ./shared_figs/temp*.png

# idealized plots
convert ./paper/fig06.png -gravity East -crop 25x20%+0+0 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -gravity East -crop 25x20%+0+0 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/dbf_ideal.png
rm -rf ./shared_figs/temp*.png
convert ./paper/fig06.png -gravity West -crop 25x20%+0+0 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -gravity West -crop 25x20%+0+0 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/dbf.png
rm -rf ./shared_figs/temp*.png


# idealized plots
convert ./paper/fig06.png -crop 640x480+0+480 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -crop 640x480+0+480 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/csh.png
rm -rf ./shared_figs/temp*.png
convert ./paper/fig06.png -gravity NorthEast -crop 640x480+0+480 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -gravity NorthEast -crop 640x480+0+480 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/csh_ideal.png
rm -rf ./shared_figs/temp*.png


convert ./paper/fig06.png -crop 640x480+0+1440 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -crop 640x480+0+1440 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/enf.png
rm -rf ./shared_figs/temp*.png
convert ./paper/fig06.png -gravity NorthEast -crop 640x480+0+1440 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -gravity NorthEast -crop 640x480+0+1440 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/enf_ideal.png
rm -rf ./shared_figs/temp*.png

convert ./paper/fig06.png -crop 640x480+0+1920 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -crop 640x480+0+1920 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/gra.png
rm -rf ./shared_figs/temp*.png
convert ./paper/fig06.png -gravity NorthEast -crop 640x480+0+1920 ./shared_figs/temp_1.png
convert ./paper/fig06b.png -gravity NorthEast -crop 640x480+0+1920 ./shared_figs/temp_2.png
convert +append ./shared_figs/temp*.png ./shared_figs/gra_ideal.png
rm -rf ./shared_figs/temp*.png

wget https://evolution.berkeley.edu/admin/media/3/4196_evo_resources_resource_image_369_original.gif

mv 4196*.gif ./shared_figs/stomata.gif
convert ./shared_figs/stomata.gif ./shared_figs/stomata.png
