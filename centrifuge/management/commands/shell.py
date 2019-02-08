from django.core.management.commands.shell import Command as ShellCommand
import os

# Command to add Ipython to the shell_plus extension
class Command(ShellCommand):
    shells = ShellCommand.shells.append('ptipython')

    def ptpython(self, arg):
        from ptpython.repl import embed
        embed(globals(), locals(), vi_mode=False)

