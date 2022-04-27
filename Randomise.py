# randomisation script
# RScript by Jakub Scaber
# Year: 2022
# Randomises all files in a folder +/- subfolders and
# outputs a dictionary, with options for filtering by filename and extension


import os
import random
import sys
import csv
import shutil
import argparse

def batch_open_images(path, file_type=None, name_filter=None, recursive=False):
    '''Open all files in the given folder.
    :param path: The path from were to open the images. String and java.io.File are allowed.
    :param file_type: Only accept files with the given extension (default: None).
    :param name_filter: Only accept files that contain the given string (default: None).
    :param recursive: Process directories recursively (default: False).
    '''

    def check_type(string):
        '''This function is used to check the file type.
        It is possible to use a single string or a list/tuple of strings as filter.
        This function can access the variables of the surrounding function.
        :param string: The filename to perform the check on.
        '''
        if file_type:
            # The first branch is used if file_type is a list or a tuple.
            if isinstance(file_type, (list, tuple)):
                for file_type_ in file_type:
                    if string.endswith(file_type_):
                        # Exit the function with True.
                        return True
                    else:
                        continue
            # The second branch is used if file_type is a string.
            elif isinstance(file_type, string):
                if string.endswith(file_type):
                    return True
                else:
                    return False
            return False
        else:
            return True

    def check_filter(string):
        '''This function is used to check for a given filter.
        It is possible to use a single string or a list/tuple of strings as filter.
        This function can access the variables of the surrounding function.
        :param string: The filename to perform the filtering on.
        '''
        if name_filter:
            # The first branch is used if name_filter is a list or a tuple.
            if isinstance(name_filter, (list, tuple)):
                for name_filter_ in name_filter:
                    if name_filter_ in string:
                        return True
                    else:
                        continue
            # The second branch is used if name_filter is a string.
            elif isinstance(name_filter, string):
                if name_filter in string:
                    return True
                else:
                    return False
            return False
        else:
            return True

    # We collect all files to open in a list.
    path_to_images = []
    # Replacing some abbreviations (e.g. $HOME on Linux).
    path = os.path.expanduser(path)
    path = os.path.expandvars(path)
    # If we don't want a recursive search, we can use os.listdir().
    if not recursive:
        for file_name in os.listdir(path):
            full_path = os.path.join(path, file_name)
            if os.path.isfile(full_path):
                if check_type(file_name):
                    if check_filter(file_name):
                        path_to_images.append(full_path)
    # For a recursive search os.walk() is used.
    else:
        for directory, dir_names, file_names in os.walk(path):
            for file_name in file_names:
                full_path = os.path.join(directory, file_name)
                if check_type(file_name):
                    if check_filter(file_name):
                        path_to_images.append(full_path)

    return path_to_images

def split_string(input_string):
    '''Split a string to a list and strip it
    :param input_string: A string that contains semicolons as separators.
    '''
    string_splitted = input_string.split(';')
    strings_striped = [string.strip() for string in string_splitted]
    return strings_striped

if __name__ == '__main__':
    #Setup the argument parser class
    parser = argparse.ArgumentParser(prog='Randomisation Script',
                                     description='''\
            Randomises all files in a folder +/- subfolders and
            outputs a dictionary,
            written by Jakub Scaber
             ''')
    parser.add_argument('--dir', action='store', help='Directory', required=True)
    parser.add_argument('--types', action='store', help='Comma separated file extensions [e.g. jpg]', default="")
    parser.add_argument('--filters', action='store', help='Comma separated strings [e.g. Mutant,Control1]', default="")
    parser.add_argument('--recursive', action='store', help='Scan subdirectories [True/False]', default=True)
    parser.add_argument('--mode', action='store', help='Whether to ', default=True)


    #Run the argument parser
    args = parser.parse_args()
    
    #Extract our value or default
    import_dir = args.dir
    file_types = args.types
    filters = args.filters
    do_recursive = args.recursive

	#Generates a list of files to open
    image_paths = batch_open_images(import_dir,
                               split_string(file_types),
                               split_string(filters),
                               do_recursive
                              )
	#Cycle through the images
    dictionary = dict()
    for image in image_paths:
        random_string = ""
        for _ in range(20):
            random_string += chr(random.randint(97, 97 + 26 - 1))
        dictionary[random_string] = image
        shutil.copyfile(image, random_string + ".czi")
    
    with open('output.csv','w') as out:
        writer = csv.writer(out)
        writer.writerows(dictionary.items())

        
        
