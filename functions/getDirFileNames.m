function filenames = getDirFileNames(directory)
%Get a cell array with file names in a specific directory


listing = dir(directory);
filenames = {listing.name};
filenames(strcmp(filenames,'.') | strcmp(filenames,'..') | strcmp(filenames,'.DS_Store')) = [];


end

