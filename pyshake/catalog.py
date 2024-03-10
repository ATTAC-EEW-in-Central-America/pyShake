from pandas import read_csv
from obspy.clients.fdsn import Client

def get_eventindex(url,method='pandas',**opt):
    """Rapid event index list from FDSN web service.

    Uses the event FDSN web service, returns the list of events.

    .. code:: python

        from pyshake.catalog import get_eventindex
        get_eventindex('http://eida.ethz.ch', limit=100)

    :param url: The fdsnws base url (e.g. 'http://eida.ethz.ch').
    :type string: :py:class:`string`.

    Any additional keyword arguments will be passed to the webservice as
        additional arguments. Passing any non-default parameters that the webservice does not support will raise an error.

    :return: Pandas dataframe.
    :rtype: ..
    """
    if method == 'pandas':

        url += '/fdsnws/event/1/query?format=text&'
        url += '&'.join(['%s=%s'%(k,opt[k])for k in opt])
        print(url)

        return read_csv(url, delimiter="|")
    
    elif method == 'obspy':

        return Client("IRIS").get_events(format='text',**opt)
    
    