import yaml

filename = "/home/yujq/users/caijie/PostProcessing/FarRadio_cmake/tests/screen.yaml"

with open(filename, "r") as f:
    data = yaml.load(f, Loader=yaml.FullLoader)
    print(data)