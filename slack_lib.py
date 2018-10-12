from slackclient import SlackClient
from configparser import ConfigParser

def print_to_slack(text, pic=None):
    c_parser = ConfigParser()
    try:
        c_parser.read('slack.cfg')
    except Exception:
        raise Exception('Config File is missing!!!!')
    login = dict()
    for key in c_parser['Basic'].keys():
        login[key] = c_parser['Basic'][key]

    args = {'channel': 'test_channel',
            'text': text}

    sc = SlackClient(**login)
    x = sc.api_call("chat.postMessage", **args)
    if pic is not None:
        y = sc.api_call('files.upload', file=open(pic, 'rb'), channels=['test_channel'])
    return x
