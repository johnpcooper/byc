from byc import plasmids, files

filepath = files.select_file("Choose the nanodrop data .csv you want to modify")

df = plasmids.extend_nanodrop_df(filepath)