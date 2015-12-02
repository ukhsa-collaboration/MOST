from subprocess import CalledProcessError
from copy import deepcopy
import uk.gov.phe.phe_exceptions
from uk.gov.phe import phe_exceptions


def validate_paths_and_fq( args, glob_for, dir_name=None ):
    '''
    Generic validation of paths to fastq files for tools using the arg_parse 
    boilerplate common to qa_and_trim, virus_variants

    Args:
        args, Namespace object : the return from parser.parse_args()
        glob_for, string : param to glob for handle to fastq[.gz] files
        dir_name, string : terminal dir appended to output_dir path

    Return:
        False, Boolean : if any arg fails to validate, otherwise
        valid_args, hash : hash of valid args
            <k> output_dir : string, output_dir path, terminates with /dir_name if not None
            <k> fq_files :  list, 2 fully-specified fastq[.gz] files

    MGGoulden 20130920
    bit of a monster func .... 123 lines 
    '''

    import os 
    import glob
    import logging
    import log_writer

    logger = logging.getLogger( 'stdout_stderr_logger' )

    invalid_args = None
    valid_args={}
    fq_files = []

    # validate input_dir && fastq[.gz] files    
    # no input_dir; applies solely to 2fastq entry point without opt --input_dir;  
    # check for fully-specified fastq paths
    if args.input_dir is None:
        pass_a_dir = '\nan input_dir may need to be passed using the -i/--input_dir flag'
        if not os.path.isfile( args.fastq_1):
            invalid_args = True
            log_txt = 'the fastq_1 argument passed ('+args.fastq_1+') does not unambiguously identify a valid file ... AMEND THIS'+pass_a_dir
            log_writer.write_log( logger, log_txt, 'info' )
        else:
            fq_files.append( args.fastq_1 )

        if not os.path.isfile( args.fastq_2 ):
            invalid_args = True
            log_txt = 'the fastq_2 argument passed ('+args.fastq_2+') does not unambiguously identify a valid file ... AMEND THIS'+pass_a_dir
            log_writer.write_log( logger, log_txt, 'info' )
        else:
            fq_files.append( args.fastq_2 )
            valid_args['fq_files'] = fq_files

    # 2fastq with opt --input_dir
    elif args.subparser_name == '2fastq':
        # does it combine with passed fq to point to file?
        if os.path.isfile( args.fastq_1 ):
            potential_fq1 = args.input_dir + '/' + os.path.split(args.fastq_1)[1]
        else:
            potential_fq1 = args.input_dir + '/' + args.fastq_1
        if not os.path.isfile( potential_fq1 ):
            invalid_args = True
            log_txt = 'the passed arguments ('+potential_fq1+') do not identify a valid fastq file ... AMEND THIS' 
            log_writer.write_log( logger, log_txt, 'info' )
        else:
            fq_files.append(potential_fq1)

        if os.path.isfile( args.fastq_2 ):
            potential_fq2 = args.input_dir + '/' + os.path.split(args.fastq_2)[1]
        else:
            potential_fq2 = args.input_dir + '/' + args.fastq_2
        if not os.path.isfile( potential_fq2 ):
            invalid_args = True
            log_txt = 'the passed arguments ('+potential_fq2+') do not identify a valid fastq file ... AMEND THIS' 
            log_writer.write_log( logger, log_txt, 'info' )
        else:
            fq_files.append(potential_fq2)
            valid_args['fq_files'] = fq_files

    # workflow & 2fastqdir have positional (mandatory) input_dir args
    else:
        fq_files = sorted( glob.glob( args.input_dir +'/'+ glob_for ) )
        if not len(fq_files) == 2:
            invalid_args = True
            log_writer.write_log( logger, 'number of fastq[.gz] files in input directory passed is not 2  ...     AMEND THIS', 'info' )
        else:
            valid_args['fq_files'] = fq_files
            log_writer.write_log( logger, 'found 2 fastq[.gz] files in the input directory passed', 'info' )


    # get output_dir
    # any of workflow, 2fastq & 2fastdir with --outdir_dir passed
    if not args.output_dir is None: 
        if not os.path.isdir( args.output_dir ):
            invalid_args = True
            log_txt = 'output_dir passed ('+args.output_dir+') is invalid ... AMEND THIS '
            log_writer.write_log( logger, log_txt, 'info' )
        else:
            valid_args['output_dir'] = args.output_dir

    # otherwise get it from elsewhere
    else: 
        # workflow & 2fastqdir defaults; 2fastq with --input_dir
        if not args.input_dir is None:
            if not os.path.isdir( args.input_dir ):
                invalid_args = True
                log_txt = 'input_dir passed ('+args.input_dir+') is invalid ... AMEND THIS '
                log_writer.write_log( logger, log_txt, 'info' )
            else:
                valid_args['output_dir'] = args.input_dir

        # 2fastq without --input_dir; get it from fq paths
        else:
            potential_out = os.path.split(args.fastq_1)[0]
            if not os.path.isdir( potential_out ):
                invalid_args = True
                log_txt = 'passed parameter ('+args.fastq_1+') does not include a valid path ... AMEND THIS;'
                log_writer.write_log( logger, log_txt, 'info' )

            potential_out = os.path.split(args.fastq_2)[0]
            if not os.path.isdir( potential_out ):
                invalid_args = True
                log_txt = 'passed parameter ('+args.fastq_2+') does not include a valid path ... AMEND THIS;'
                log_writer.write_log( logger, log_txt, 'info' )

            else:
                valid_args['output_dir'] = potential_out
                
    # return     
    if invalid_args:
        return False
    else:
        if dir_name:
            valid_args['output_dir'] = valid_args['output_dir']  +'/'+ dir_name
        return valid_args




def get_from_config( CONFIG_FILE, get_items=None ):
    '''  
    Returns data for the calling function from CONFIG_FILE. 

    Args
        CONFIG_FILE, string : environment variable loaded using module, pointing to config file
        get_items, list : optional, list of items to return from config file

    Return
        data-structure : default return is the data structure represented by the config file for the calling function
        items, unpacked tuple : if get_items is specified, unpacked variables specified by the config file for the calling function 

    The config_file.yml must represent a dictionary, with the primary key == name_of_func_using_value
    <k> = function_using_FOO
        <k> = FOO
    `       <v> = value_of_foo
    'function_using_FOO' is retrieved from the call stack and used selectively to return data from the config_file.yml
    
    N.B. the calling code for a single item needs a comma to unpack the return from tuple
    return_item, = get_from_config( config_file, get_items=['return_item'] )

    MGGoulden 20130906
    amended 20130919
    '''

    import os
    import yaml
    import logging
    import inspect 
    # from common_modules
    import log_writer # setup_logger, write_log, error_header, info_header, log_process    


    # get a pointer to the logger
    logger = logging.getLogger( 'stdout_stderr_logger' ) 
    log_txt = 'IN get_from_config( CONFIG_FILE = '+CONFIG_FILE+', get_items = ' +str(get_items)+ ')'
    log_writer.info_header( logger, log_txt )

    # sanity check
    if not CONFIG_FILE in os.environ:
        log_txt = 'The environment variable ('+CONFIG_FILE+') is not available; module load may be required; quitting ... '
        log_writer.error_header( logger, log_txt )
    elif not os.path.exists( os.environ[ CONFIG_FILE ] ):
        log_txt = 'The environment variable ('+CONFIG_FILE+') points to a file ('+os.environ[CONFIG_FILE]+') which does not exist; module amendment required; quitting ... '
        log_writer.error_header( logger, log_txt )

    # who calls? - identify the calling func from stack
    caller =  inspect.stack()[1][3]
    if caller == 'try_and_except':
        caller =  inspect.stack()[2][3]

    read_me = os.environ[CONFIG_FILE]    
    with open( read_me, 'r') as f:
        CONFIG = yaml.load(f)
        if not get_items:
            return CONFIG[caller]
        else:
            return [ CONFIG[ caller ][x] for x in get_items ]




def try_and_except(error_filepath, function, *parameters, **named_parameters):
    """
    This wraps a function in try and except clause. If an error is caught this will trigger
    1) reporting of the error to stdout
    2) writing of the error into an error file
    3) a sys.exit with code 1

    Parameters
    ----------
    error_file : String
        path to log file to capture error in
        (a useful default = logger.handlers[0].baseFilename, which returns 
            the stderr.log FileHandler added first by log_writer.setup_logger)
    function: Function
        the function
    parameters : all non-named parameters for the function
    named_parameters : all named paramaters for the function

    Notes
    -----
    Returns the returns from the function called

    Examples
    --------

    assuming a function

    def my_func(a,b,c = None)
        ......
        return d

    This function would be called as follows:

    return_val = try_and_except("stderr.log", my_func, 1, 2, c = 3)
    """
    import sys, traceback, subprocess, logging
    try:
        return function(*parameters, **named_parameters)
    except phe_exceptions.PheException as phe_e:
        # This exception is created by the 'call_external', when exit code != 0
        
        logger = logging.getLogger("stdout_stderr_logger")
        logger.exception(function.__name__ + " has raised an exception: \n" + str(phe_e))
        
        # Exit with the returncode specified in the called process.
        # TODO: Should it be a different return code? E.g. ranged for traceback.
        sys.exit(phe_e.phe_returncode)
    except Exception as e:
        error_string = "There was an error in the function '" + function.__name__ + "'"
        error_divider = "_" * 60
        print error_string
        print error_divider
        traceback.print_exc()
        print error_divider

        error_file = open(error_filepath, "a")
        error_file.write(error_string + "\n")
        error_file.write(error_divider + "\n")
        traceback.print_exc(file = error_file)
        error_file.write(error_divider  + "\n")
        error_file.close()
        sys.exit(1)

def check_file_exists(filepath, file_description):
    """
    A function to check if a file exists.
    It will print out an error message and exit if the file is not found

    Parameters
    ----------
    filepath : String
        the path to the file to be checked
    file_description : String
        a description of the file to be checked e.g "config file"
    """
    import os, sys
    if not os.path.exists(filepath):
        print("The " + file_description + " (" + filepath + ") does not exist")
        sys.exit(1)

def add_module_dir_to_path(module_dir):
    import sys
    if module_dir not in sys.path:
        sys.path.insert(1, module_dir)

def add_module_dir_relative_path_to_sys_path(relative_path):
    import os, inspect
    module_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],relative_path)))
    add_module_dir_to_path(module_dir)

def add_standard_module_dirs_to_path():
    standard_module_dirs = ["modules", "../common_modules", "../../common_modules"]
    for dir in standard_module_dirs:
        add_module_dir_relative_path_to_sys_path(dir)

def add_workflow_module_dirs_to_path(workflow, config_file_path = "/Volumes/NGS2_DataRAID/software/workflows/config_files/"):
    import yaml
    yaml_file = open(config_file_path)
    config = yaml.load(yaml_file)
    module_dir_relative_paths = config[workflow]
    for module_dir_relative_path in module_dir_relative_paths:
        add_module_dir_relative_path_to_sys_path(module_dir_relative_path)
        
def call_external(cmd, logger, raise_exception=False):
    '''
    Calls specified external command, waits for the process to finish and returns
    whatever the child process has returned. If 'raise_exception' is True
    a CalledProcessError exception is raised, which includes the return code,
    command and the output (stdout + stderr) of the external process. The call can be 
    wrapped inside the try_and_except function for handling the possible raised 
    exceptions (should be used with 'raise_exception'=True). Otherwise, the
    exception can be handled manually.
    
    @param cmd: Command to be ran with all appropriate arguments. The validity of
        the command the the arguments are not checked.
    @type cmd: arr.
    @param logger: Logger to be used for logging the output from the process.
    @type logger: logger.
    @param raise_exception: Specifies whether a CalledProcessError should be raised
        when returncode is not 0 (default False).
    @type raise_exception: bool.
    
    @return: Returns the returncode from the external process.
    
    @raise CalledProcessError: If returncode is not 0 and 'raise_exception' is
        set to True, then CalledProcessError will be raised.
    '''
    import subprocess
    import log_writer
    
    # For now (v 2.7) can't use subprocess.call, subprocess.check_all because PIPE
    #    is not correctly reading the method. For now, use Popen and wait().
    process = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE) 
    process.wait()
    
    # Log the outputs of the external programm.
    # FIXME: 'log_error_to' will always go to 'error' if exit code > 0 (see log_process) 
    process_out, process_err = log_writer.log_process(logger, process, log_error_to = "info")
    
    # Use stdout as output, unless err is not empty. In which case append it.
    # This may pollute the error log with stdout form the process.
    out = deepcopy(process_out)
    if process_err:
        out = out.join(["********ERROR********\n", process_err])
    
    if raise_exception and process.returncode != 0:
        raise phe_exceptions.PheExternalError("External script has returned non-zero exit code.", 
                                        subprocess.CalledProcessError(process.returncode, cmd, out))
    else:
        retval={'proc_returncode' : process.returncode, 'proc_stdout' : process_out, 'proc_stderr' : process_err }
        return retval

def write_component_complete( output_dir ):
    '''
    Creates marker of complete module run by writing empty file 'module_complete.txt' to result/module/.
    Intended to be called as final call from module.

    Args
        output_dir, str : full path to the component output dir    

    Return
        none

    Side-effects
        writes empty file 'ComponentComplete.txt' into output_dir
    '''
    cc = '/'.join([ output_dir, 'ComponentComplete.txt' ])
    with open( cc, 'w' ) as cc_out:
        cc_out.write('')


