"""
Create single redis instance and import it around where it is needed
"""
import redis
redis_instance = redis.StrictRedis(
    host="127.0.0.1", port=6379, db=0, decode_responses=True
)
