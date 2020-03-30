import sys
import argparse
import json

import requests


def execute_command_as_string(data, host=None, port=None):
    """"""
    r = requests.post(
        "http://{}:{}/jsonrpc".format(host, port),
        data=data,
        headers={"Content-Length": str(len(data)), "Content-Type": "application/json"},
    )
    try:
        json_respond = json.loads(r.text)
        return json_respond
    except Exception as err:
        # FIXME: raise
        print(err)


def send_message_port(message, ip_address, port):
    message_to_send = (
            '{"id":"1", "method":"user_message","params":{"content":"%s"}}' % message
    )
    results = ""
    try:
        results = execute_command_as_string(message_to_send, ip_address, port)
    except Exception as err:
        # FIXME: raise
        print("message send fail", err)
    return results



def main(parser, args):
    while 1:
        msg = input("Enter a message to send: ")
        send_message_port(msg, args.host, args.port)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Send messages to MinKNOW. Run `mk_manager_client --list` and choose the port "
                        "for your device."
    )
    parser.add_argument("--host", default="127.0.0.1", help="default: 127.0.0.1")
    parser.add_argument("--port", type=int, default="8007", help="defaul: 8007")
    try:
        main(parser, parser.parse_args())
    except KeyboardInterrupt:
        sys.exit(0)
