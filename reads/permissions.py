from rest_framework.permissions import BasePermission


class IsOwner(BasePermission):
    """
    Custom permission to only allow owner of an object to access it.
    """

    def has_object_permission(self, request, view, obj):
        return obj.owner == request.user
