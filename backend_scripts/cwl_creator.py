from configparser import ConfigParser
import classes
import os
import sys

config = ConfigParser()
print(os.getcwd())
os.chdir(os.path.join(os.path.dirname(__file__), '..' ))
print(os.getcwd())
config.read("../config.ini")
print(config.get("main", "database"))
test = classes.database_checker(config.get("main", "database"))
test.check_and_run(config.get("main", "root"))
