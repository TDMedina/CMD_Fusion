
import argparse
import sys


class CustomHelp(argparse.HelpFormatter):
    """Custom help formatter_class that only displays metavar once."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar
        parts = []
        if action.nargs == 0:
            parts.extend(action.option_strings)
        else:
            default = self._get_default_metavar_for_optional(action)
            args_string = self._format_args(action, default)
            for option_string in action.option_strings:
                parts.append("%s" % (option_string))
            parts[-1] += " %s " % args_string
        return ", ".join(parts)

    def _format_action(self, action):
        parts = super()._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts

def _setup_argparser():
    argparser = argparse.ArgumentParser("FusionDash")
    argparser.add_argument("-g", "--gtf", required=True)
    argparser.add_argument("-s", "--sample-dir", required=True)
    argparser.formatter_class = CustomHelp
    return argparser

def parse_dashboard_args():
    parser = _setup_argparser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    args = vars(args)
    return args
