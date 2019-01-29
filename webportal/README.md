# Random notes on Django

## local set up
##### NOTE: the app's migration folder and the sqlite database are not commited to keep database separated from the rest (see .gitignore)

To setup database migration:
```
$ python manage.py makemigrations
$ python manage.py migrate
$ python manage.py runserver
```

## identifier
see ./analysis/models.py

use python's uuid module to add unique identifier to each session


## choices
when defining models you can set choices for the CharField
