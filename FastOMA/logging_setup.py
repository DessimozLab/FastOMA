import logging


def setup_logging(verbosity: int, thread_names: bool = False):
    """
    Configure logging based on verbosity levels:
    -v   : INFO level on main package, WARNING on sub-packages
    -vv  : DEBUG level on main package, INFO on sub-packages
    -vvv : DEBUG on all packages
    """
    log_levels = {
        0: (logging.WARNING, logging.WARNING),  # Default: WARNING everywhere
        1: (logging.INFO, logging.WARNING),     # -v
        2: (logging.DEBUG, logging.INFO),       # -vv
        3: (logging.DEBUG, logging.DEBUG),      # -vvv
    }

    # Get the appropriate levels
    main_level, sub_level = log_levels.get(max(0, verbosity), (logging.DEBUG, logging.DEBUG))

    # Set main package logger
    logging.getLogger("FastOMA").setLevel(main_level)

    # Set sub-packages
    logging.getLogger("FastOMA.zoo").setLevel(sub_level)
    if thread_names:
        fmt = "%(asctime)s %(levelname)-8s [%(threadName)s] %(name)s: %(message)s"
    else:
        fmt = "%(asctime)s %(levelname)-8s %(name)s: %(message)s"

    # Configure root logger
    logging.basicConfig(
        format=fmt,
        level=min(main_level, sub_level),  # Ensure root captures the lowest level needed
    )
