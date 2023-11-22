name_db="combined_rcc_ivt_train_pqlist"
model_loc="../resource/model_ivt_only_pub"
python train_ivt_48.py $name_db $model_loc
python train_ivt_with_aug_48.py $name_db $model_loc
