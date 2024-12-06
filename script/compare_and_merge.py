import re
def compare_files(file_1, file_2, ):
    with open(file_1, 'r') as file:
        dat_1 = file.readlines()

    with open(file_2, 'r') as file:
        dat_2 = file.readlines()

    for txt in dat_2:
        if txt not in dat_1:
            print(f"adding {txt}")
            dat_1.append(txt)   

    return(dat_1)

if __name__ == "__main__":
    full_txt = compare_files("data/glider/filipa_txt/biocarbon_glider_standard_sensors_all.txt", "data/glider/filipa_txt/noc-grease409_sensors.txt")
    with open('full_list.txt', 'w') as f:
        for line in full_txt:
            f.write(line)