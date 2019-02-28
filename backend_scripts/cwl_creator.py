from configparser import SafeConfigParser
import classes

config = SafeConfigParser()
config.read("../config.ini")
print(config.get("main", "database"))
test = classes.database_checker(config.get("main", "database"))
test.check_and_run(config.get("main", "root"))
