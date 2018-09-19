from django.core.management import BaseCommand, CommandError

from web.tasks import update_flowcell_details


class Command(BaseCommand):

    help = 'Run update_flowcell_details task'

    def handle(self, *args, **options):

        try:

            print('Running update_flowcell_details task')

            update_flowcell_details()

        except Exception as e:

            raise CommandError(repr(e))
