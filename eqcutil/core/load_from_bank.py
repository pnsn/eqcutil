import logging

from eqcorrscan import Tribe
from obsplus import WaveBank, EventBank

from eqcutil.augment.catalog import apply_phase_hints, filter_picks, transfer_picks
from eqcutil.augment.template import rename_templates
from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import rich_error_message

Logger = logging.getLogger(__name__)

def generate_clustering_tribe_from_banks(
    wavebank: WaveBank,
    eventbank: EventBank,
    event_ids: list,
    pick_filt_kwargs={'enforce_single_pick': 'preferred'},
    creation_kwargs={'method':'from_client'},
    transfer_mapping={},
    rename=True) -> ClusteringTribe:

    # Compatability check for wavebank
    if not isinstance(wavebank, WaveBank):
        raise TypeError('wavebank is not type obsplus.bank.wavebank.WaveBank')
    elif len(wavebank.read_index()) == 0:
        raise ValueError('wavebank is empty')
    
    # Ensure wavebank is mapped to creation kwargs
    creation_kwargs.update({'client_id': wavebank,
                            'method': 'from_client'})
    # Compatability check for eventbank
    if not isinstance(eventbank, EventBank):
        raise TypeError('eventbank is not type obsplus.bank.eventbank.EventBank')
    elif len(eventbank.read_index()) == 0:
        raise ValueError('eventbank is empty')
    # Read indices
    wbidx = wavebank.read_index()
    ebidx = eventbank.read_index()
    
    # Compatability check for event_ids
    if hasattr(event_ids, '__iter__'):
        ebidx_sub = ebidx[ebidx.event_id.isin(event_ids)]
        if len(ebidx_sub) == 0:
            raise ValueError('no matches of event_ids in eventbank')
    else:
        raise AttributeError('event_ids must be an iterable comprising QuakeML event_id strings')
    
    # PROCESSING
    ctr = ClusteringTribe()
    for event_id in ebidx_sub.event_id:
        cat = eventbank.get_events(event_id=event_id)
        cat = apply_phase_hints(cat)
        cat = transfer_picks(cat, transfer_mapping)
        cat = filter_picks(cat, **pick_filt_kwargs)
        if len(cat) == 0:
            continue
        if len(cat[0].picks) == 0:
            continue
        creation_kwargs.update({'catalog': cat})
        try:
            tribe = Tribe().construct(**creation_kwargs)
        except Exception as e:
            Logger.warning(f'{event_id} produced {rich_error_message(e)}')
            continue
        # Ensure template stream is merged
        for template in tribe:
            template.st.merge()
        # Apply renaming if specified
        if rename:
            tribe = rename_templates(tribe)
        # Apply template mapping to channels
        # TODO: enforce rules to only allow location and component changes
        for template in tribe:
            for tr in template.st:
                if tr.id in transfer_mapping.keys():
                    new_id = transfer_mapping[tr.id]
                    new_id_parts = new_id.split('.')
                    for _e, _k in enumerate(['network','station','location','channel']):
                        if tr.stats[_k] != new_id_parts[_e]:
                            tr.stats[_k] = new_id_parts[_e]
        ctr += tribe

    return ctr

