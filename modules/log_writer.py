import logging

def setup_logger(info_file = "stdout.log", error_file = "stderr.log", logger_name = 'stdout_stderr_logger' ):
    """
    A function to set up a logger for writing out log files

    Parameters
    ----------
    info_file : String
        Path to info level log file (default: stdout.log)
    error_file : String
        Path to error level log file (default: stderr.log)

    Returns
    -------
    logger : A instance of the Logger class
    """

    class LogLevelFilter(object):
        def __init__(self, level):
            self.__level = level

        def filter(self, logRecord):
            return logRecord.levelno <= self.__level

    logger = logging.getLogger( logger_name )
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s\n%(message)s')

    handler_stderr = logging.FileHandler(error_file)
    handler_stderr.setLevel(logging.ERROR)
    handler_stderr.setFormatter(formatter)
    handler_stderr.addFilter(LogLevelFilter(logging.ERROR))
    logger.addHandler(handler_stderr)

    handler_stdout = logging.FileHandler(info_file)
    handler_stdout.setLevel(logging.INFO)
    handler_stdout.setFormatter(formatter)
    handler_stdout.addFilter(LogLevelFilter(logging.INFO))
    logger.addHandler(handler_stdout)

    return logger

def write_log(logger, log_text,log_level):
    """
    Writes text to a logger at a particular log level

    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    log_text : String
        The text to be written to the log file
    log_level : String
        The level of logging to which the text should be written (either 'info' or 'error')

    """

    if log_level == "error":
        logger.error(log_text)
    elif log_level == "info":
        logger.info(log_text)

def log_process(logger, process,log_info_to = "info", log_error_to = "error", limit_logging = 0):
    """
    A function to log the output of a subprocess.Popen call

    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    process : process pipe
        A process created by subprocess.Popen
    log_info_to: String
        The level at which to log info level logs into (default 'info')
    log_error_to: String
        The level at which to log error level logs into (default 'error')
    limit_logging: Integer
        Limit logging to either stdout (1) or stderr (2), log both if 0 (default 0)

    Returns
    -------
    stdout : String
        The stdout from the process
    stderr : String
        The stderr from the process
    """
    #if the process gives a exit status greater than 1, i.e. an genuine error, then log_error_to 'error'.
    if process.returncode > 0:
        log_error_to = "error"

    stdout = ""
    stderr = ""

    if limit_logging != 2:
        if not process.stdout == None:
            stdout = process.stdout.read()
        if len(stdout) > 0:
            write_log(logger, stdout, log_info_to)
    if limit_logging != 1:
        if not process.stderr == None:
            stderr = process.stderr.read()
        if len(stderr) > 0:
            write_log(logger, stderr, log_error_to)

    return stdout, stderr

def write_header_to_log(logger, header_text, log_level):
    """
    A utility function to write some header text bound by asterisks to the log
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    log_level : String
        The level at which to log
    """
    write_log(logger, "******** " + header_text + " ********", log_level)

def error_header(logger, header_text):
    """
    A utility function to write some header text bound by asterisks to the error log_text
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    """
    write_header_to_log(logger, header_text, "error")

def info_header(logger, header_text):
    """
    A utility function to write some header text bound by asterisks to the info log_text
    Parameters
    ----------
    logger : Logger
        An instance of the Logger class
    header_text : String
        The header text to embed within asterisks
    """
    write_header_to_log(logger, header_text, "info")

    
def get_logger_path( args, dir_name='logs' ):
    ''' 
    Returns a path for logger output according to params in the return from
    parser.parse_args(). By default returns output_dir/logs, else returns
    input_dir/logs, else fastqfile_dir/logs

    Args:
        args, Namespace object : the return from parser.parse_args() 
        dir_name, string : terminal dir appended to derived path

    Returns:
        output_dir, string : fully-specified output path terminating
        with /dir_name to which logs can be written

    Side effect:
        creates path if necessary   

    MGGoulden 20130709
    amended 20130927 & moved into log_writer module
    amended 20140109 to make more generic

    '''
    
    # dependencies  
    import os, sys

    log_path = None
    # if an output_dir has been passed check it & use it
    if not args.output_dir is None:
        if not os.path.isdir( args.output_dir ):
            print 'the output_dir passed ('+ args.output_dir +') does not exist'
            print 'making dir: '+ args.output_dir
            os.mkdir( args.output_dir )
        log_path = args.output_dir

    # otherwise use any input_dir
    elif not args.input_dir is None:
        if not os.path.isdir( args.input_dir ):
            print 'the input_dir passed ('+ args.input_dir +') is not valid'
        else:
            log_path = args.input_dir

    # with neither input_dir nor output_dir passed, extract a path from fastq files passed
    elif args.fastq_1 and args.fastq_r2:
        fqs = [ args.fastq_1, args.fastq_r2 ]
        fq_path = set([ os.path.split(fq)[0] for fq in fqs ])
        if not len(fq_path) == 1:
            print 'fastq files passed ('+ str(fqs) +') have different paths'
        else:
            (fq_path,) = fq_path
            if not os.path.isdir( fq_path ):
                print ''.join(['fastq files passed have invalid path: ', fq_path ])
            else:
                log_path = fq_path

    assert log_path is not None,'a valid path for logger could not be derived from the passed args'

    full_log_path = ''.join( [ log_path, '/', dir_name ] )
    if not os.path.isdir( full_log_path ):
        os.mkdir( full_log_path )
    return full_log_path
