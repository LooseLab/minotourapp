from rest_framework import serializers
from tabs.models import RunID, RunTabs
# from reads.models import MinIONRun


# class RunTabsSerializer(serializers.ModelSerializer):
class RunIDSerializer(serializers.ModelSerializer):

    class Meta:
        model = RunID
        fields = ('run_id_id', 'tab_id_id',)
