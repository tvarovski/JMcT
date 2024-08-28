import functools
import logging
from typing import Callable

logger = logging.getLogger('__main__.' + __name__)

def time_elapsed(func: Callable) -> Callable:
    """
    Decorator that logs the time taken by a function to execute.

    Parameters:
    func (Callable): The function to be timed.

    Returns:
    Callable: The wrapped function which will log its execution time.
    """
    import time
    # decorator that times the elapsed time of a function
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """
        Wrapper function that calculates and logs the execution time of the decorated function.

        Parameters:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

        Returns:
        The result of the decorated function.
        """
        #start the timer
        start = time.time()
        #call the function
        result = func(*args, **kwargs)
        #stop the timer
        end = time.time()
        #calculate the elapsed time
        elapsed = end - start
        #print the elapsed time rounded to 2 decimal places
        logging.info(f"\n{func.__name__} took {elapsed:.2f} seconds to run...")
        return result
    
    return wrapper