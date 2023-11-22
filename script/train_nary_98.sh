name_db="/mnt/share/bhaskar/rcc_data_extend/6k/species"
model_loc="../resource/model_comb_ivtmap_6k_pub/"
python train_nary.py $name_db $model_loc
python train_nary_with_aug.py $name_db $model_loc

