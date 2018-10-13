from slackclient import SlackClient
from slack_config import login
def print_to_slack(text, pic=None):
    args = {'channel': 'test_channel',
            'text': text}

    sc = SlackClient(**login)
    x = sc.api_call("chat.postMessage", **args)
    if pic is not None:
        y = sc.api_call('files.upload', file=open(pic, 'rb'), channels=['test_channel'])
    return x
