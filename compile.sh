if [ "$#" -ne 1 ]; then
  echo "Usage: sh $0 (0: CPU, 1: GPU)" >&2
  exit 1
fi

make IFT_GPU=$1 preproc bag_of_feature_points create_layer_model merge_layer_models encode_layer encode_merged_layer decode_layer detection delineation

