name_db="/mnt/share/bhaskar/rcc_data_extend/6k/all_species_rcc"
model_loc="../resource/model_comb_ivtmap_rcc_6k_pub/"
python train_nary.py $name_db $model_loc
python train_nary_with_aug.py $name_db $model_loc
