import logging
import os
from configparser import ConfigParser

import classes

if __name__ == "__main__":
	# change current working directory to project root
	root = os.path.realpath(os.path.join(os.path.dirname(__file__), '../..' ))
	os.chdir(root)

	# read config file
	config = ConfigParser()
	config.read("config.ini")
	database = config.get("main", "database")

	test = classes.database_checker(database)
	test.check_and_run(root)
