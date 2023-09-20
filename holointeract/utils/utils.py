import os.path


def create_new_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        # print(f'{dir_path} directory created')
    else:
        print(f'{dir_path} directory already exists')
