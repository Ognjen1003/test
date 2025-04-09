from loguru import logger
from typing import Optional

class LoggerSingleton:
    _instance: Optional['LoggerSingleton'] = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(LoggerSingleton, cls).__new__(cls)
            cls._configure_logger()
        return cls._instance

    @staticmethod
    def _configure_logger():
        logger.add("Logs/app.log", rotation="1 MB", retention="30 days", level="INFO")
        print('Logger configured level: app.log, INFO, retention 30 days.')


    @staticmethod
    def log_info(message: str):
        logger.info(message)
        print(f'logger info: {message}')

    @staticmethod
    def log_info(message: str):
        logger.info(message)
        print(f'logger info: {message}')
    
    @staticmethod
    def log_warning(message: str):
        logger.warning(message)
        print(f'logger warning: {message}')
   
    @staticmethod
    def log_error(message: str):
        logger.error(message)
        print(f'logger error: {message}')
    
    @staticmethod
    def log_error(message: str, ip_address: str ):
        logger.error(f"{message} - IP: {ip_address}")
        print(f"{message} - IP: {ip_address}")