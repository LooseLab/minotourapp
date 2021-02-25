import gzip
import os
from io import BytesIO

from django.core.exceptions import ImproperlyConfigured


def get_env_variable(var_name):
    """ Get the environment variable or return exception """
    try:
        return os.environ[var_name]
    except KeyError:
        error_msg = "Set the %s env variable" % var_name
        raise ImproperlyConfigured(error_msg)


class GZIPBodyDecompressMiddleware:
    """
    Decompress any gzipped bytes objects being posted or put by minFQ
    """
    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.

    def __call__(self, request):
        # Code to be executed for each request before
        # the view (and later middleware) are called.
        encoding = request.headers.get("Content-Encoding", "utf-8")

        if encoding == "gzip":
            request._body = gzip.decompress(request.body)
            request._stream = BytesIO(request.body)

        response = self.get_response(request)

        # Code to be executed for each request/response after
        # the view is called.

        return response