FROM redis:latest
RUN mkdir -p /redis_files
COPY redis.conf  /usr/local/etc/redis/redis.conf
RUN chmod -R 755  /usr/local/etc/redis/redis.conf
RUN chmod -R 766  /redis_files
CMD [ "redis-server", "/usr/local/etc/redis/redis.conf" ]