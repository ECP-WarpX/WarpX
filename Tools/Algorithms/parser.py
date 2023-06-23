def parse_input(input_file):
    """
    Parse WarpX input file.

    Parameters
    ----------
    input_file : string
        Path to input file.

    Returns
    -------
    input_dict : dictionary
        Dictionary storing WarpX input parameters
        (parameter's name stored as key, parameter's value stored as value).
    """
    input_dict = dict()
    for line in open(input_file):
        sline = line.split('=')
        # skip lines that are commented out or blank
        skip_line = sline[0].startswith('#') or sline[0].startswith('\n')
        if not skip_line:
            key = sline[0].strip()
            val = sline[1].split()
            # The value corresponding to a given key of input_dict is a list
            # of strings, from which we remove any leftover comments
            for w in val:
                if w.startswith('#'):
                    val_index = val.index(w)
                    val = val[:val_index]
            input_dict[key] = val
    return input_dict
