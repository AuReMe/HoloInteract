import os.path


def rename_id(species, path, gmove=False):
    gbk_path = os.path.join(path, species+'.gbk')
    new_gbk_dir = os.path.join(path, 'gbk_renamed')
    new_gbk = os.path.join(new_gbk_dir, species+'.gbk')
    dict_dir = os.path.join(path, 'dict')
    dict_assoc_id = os.path.join(dict_dir, species+'_dict.tsv')
    if not os.path.exists(new_gbk_dir):
        os.makedirs(new_gbk_dir)
    if not os.path.exists(dict_dir):
        os.makedirs(dict_dir)
    with open(gbk_path, 'r') as igbk, open(new_gbk, 'w') as ogbk, open(dict_assoc_id, 'w') as odic:
        for l in igbk:
            if '/locus_tag=' in l:
                id = l.split('/locus_tag=')[1]
                if gmove:
                    new_id = id[0] + id[6:]
                else:
                    new_id = id
                invalid_characters = ['-', '|', '/', '(', ')', '\'', '=', '#', '*',
                                      '.', ':', '!', '+', '[', ']', ',', " "]
                for ch in invalid_characters:
                    new_id = new_id.replace(ch, '_')
                l = l.split('/locus_tag=')[0] + '/locus_tag=' + new_id
                if id != new_id:
                    if gmove and id[1:5] == "prot":
                        odic.write(new_id[:-1] + "\t" + id)
                    elif not gmove:
                        odic.write(new_id[:-1] + "\t" + id)
            ogbk.write(l)
    if os.path.getsize(dict_assoc_id) == 0:
        os.remove(dict_assoc_id)

    return new_gbk


"""path = 'e2g/public'
species = 'Sargassum-fusiforme'
rename_id(species, path)"""
# for species in os.listdir(os.path.join(path, 'gbk')):
#     species = species[:-4]
#     rename_id(species, path)
