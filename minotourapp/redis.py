"""
Create single redis instance and import it around where it is needed
"""
import redis

from minotourapp.utils import get_env_variable

redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=1, decode_responses=True
) if "localhost" in get_env_variable("MT_DJANGO_REDIS_URL") else redis.StrictRedis(
        unix_socket_path=get_env_variable("MT_DJANGO_REDIS_URL").split("//")[-1], decode_responses=True
)

