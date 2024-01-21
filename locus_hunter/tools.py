import os


def get_temp_path(prefix: str = 'temp', suffix: str = '') -> str:
    i = 0
    while True:
        path = f'{prefix}_{i:06d}{suffix}'
        if not os.path.exists(path):
            break
        i += 1
    return path
