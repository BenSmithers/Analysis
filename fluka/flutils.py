import os # check for file existence 

"""
Ben Smithers
benjamin.smithesr@mavs.uta.edu

A few tools I wrote to work with FLUKA datafiles 
"""

# Used to skip over the headers
line_width = 10

def parse_line( what ):
    """
    Parses a line in resnuclei output file and returns a list of integers 

    Raises a TypeError if 'what' is not a string
    Raises a ValueError if it fails to parse the line 
    """
    if not isinstance(what, str):
        raise TypeError("Expected {}, got {}".format(str, type(what)))

    divided = what.split(" ")
    divided = list(filter( lambda el: el!='', divided))
    parsed = [ float( part ) for part in divided ]
    # hacky way of avoiding the header lines with only numbers 
    if len(parsed)!=line_width:
        raise ValueError()

    return( parsed )

def load_datafile( which ):
    """
    Parses the text file at 'which', assumed one of those stupid FLUKA 1D streams reformatted as a 2D one

    Reads the stuff right out and converts to floats

    Does no reshaping! 
    """
    if not os.path.isfile( which ):
        raise IOError("Couldn't find the file at {}".format(which ))

    # go line by line and parse the data
    data_array = []
    with open(which, 'rt') as obj:

        line = obj.readline()
        while line:
            line = obj.readline()
            try:
                parsed = parse_line(line)
                data_array.append( parsed )
            except ValueError:
                # exception on header files. Squash and move on
                pass

    if data_array==[]:
        raise IOError("No Data parsed?")

    return(data_array)

def list_shape( source, as_str=True ):
    """
    A bit like an analog to the numpy shape method but for regular python lists. 

    Accepts a *2D* list, `source`, and returns either a string "(MxN)" or a length-2 list representing its length
    """
    if not isinstance(source, list):
        raise TypeError("Expected {}, got {}".format(list, type(source)))
    if not isinstance(source[0], list):
        raise TypeError("Expected list contents of type {}, got {}".format(list, type(source[0])))
    
    if as_str:
        return( "{}x{}".format(len(source), len(source[0])) )
    else:
        return([ len(source), len(source[0])])


def reshape( source, new_shape):
    """
    Reshapes a FLUKA 1D datastream, which was reshaped into a 2D one, into the desired format

    source - the 2D formatted fluka table from the formatted output file
    new_shape  - a length 2 list representing the desired shape 
    """
    if not isinstance(new_shape, list):
        raise TypeError("Shape should be a {}, got {}".format(list, type(new_shape)))
    if not len(new_shape)==2:
        raise ValueError("Shape should be length-2, not {}".format(len(new_shape)))
    for i in range(2):
        if not isinstance(new_shape[i], int):
            raise ValueError("Entry {} should be {}, got {}".format(i, int, type(new_shape[i])))

    old_shape = list_shape(source, False)
        
    # the datafile resnuc puts out is a 1D stream that wraps to the next row after 10 entries
    # so we need to reshape this bullshit 
    reform = [[0.0 for z in range(new_shape[1])] for nmzmk in range( new_shape[0] )]

    # we need to keep track of our place in each 2D array independently 
    row = 0
    column = 0
    ref_row = 0
    ref_column = 0

    # and make sure there are hte same number of entries in each. This is imperative! 
    assert( (len(reform)*len(reform[0])) == (len(source)*len(source[0])))

    # so now we scan across the source (source) and assign the values to the reformed dataarray 
    while row<(len(source)):
        reform[ref_row][ref_column] = source[row][column]

        ref_column+=1
        if ref_column==len(reform[ref_row]):
            ref_column = 0
            ref_row +=1

        column+=1
        if column==len(source[row]): # line_width:
            column = 0
            row+=1

    return( reform )

