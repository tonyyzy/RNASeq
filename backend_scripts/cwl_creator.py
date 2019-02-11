import programs
import classes
import os

print(os.getcwd())
ale = classes.database_reader(1,"2")
ale.extract_from_database()
print(ale.Index)
print(ale.Mapper)
print(ale.Assembler)
print(ale.Analysis)
print(ale)
ale_logic = classes.logic_builder()
ale_logic.create_workflow_logic(ale)
ale_logic.create_indexing(ale)
