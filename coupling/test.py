def InitializeCoupledSimulation(path_to_ctrl, debug):
    '''
    Function to set all required default data, e.g. well names, paths, ...
    :param path_to_ctrl: path to input file containing control data
    :param type: str
    :returns: no return value
    '''
    title_loc = ''
    dir_loc = ''

    if path_to_ctrl[-1] == '\\':
        del path_to_ctrl[-1]
    
    idxs = [i for i,key in enumerate(path_to_ctrl) if key=='\\']

    if not path_to_ctrl[-10:] == '.main_ctrl':
        title_loc = path_to_ctrl[idxs[-1] + 1:]
        dir_loc = path_to_ctrl[:-len(title_loc)]
        path_to_ctrl += '\\.main_ctrl'
    else:
        print(idxs)
        pos_diff = len(path_to_ctrl) - idxs[-1] + 1
        print(len(path_to_ctrl), pos_diff)
        title_loc = path_to_ctrl[idxs[-1] + 1:-(pos_diff - 10)]
        dir_loc = path_to_ctrl[:idxs[-1] + 1]


    print('Reading inputile: ', title_loc)
    print('in working directory: ', dir_loc)

