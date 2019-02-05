docker run --name minotour-redis -v /Users/roberto/projects/minotour/data/redis_data:/data -d redis redis-server --appendonly yes

python manage.py runserver 0.0.0.0:8000
