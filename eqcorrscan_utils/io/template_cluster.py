"""
Supporting scripts to help prepare and save inputs and outputs
of cross-correlation based clustering of :class:`~eqcorrscan.Template`
objects
"""

import os, logging
from pathlib import Path


Logger = logging.getLogger(__name__)

# def make_template_list(tribe, use_template_names=False):
#     if use_template_names:
#         return [(_t.st, _t.name) for _t in tribe]
#     else:
#         return [(_t.st, _e) for _e, _t in enumerate(tribe)]
    
def save_cluster_results(tribe, groups, savedir='cluster_results', save_subtribes=False):
    """Save the results of a :meth:`~eqcorrscan.util.clustering.cluster` call that had
    the distance matrix saved to the current working direcory (dist_mat.npy) and used
    the decorator :meth:`~eqcorrscan_utils.util.decorators.save_args` to save the
    method's input arguments to csv (cluster_params.csv).

    This method creates a save directory, moves these two files into it, and then
    uses the **groups** output from :meth:`~eqcorrscan.util.clustering.cluster`
    and the originating :class:`~eqcorrscan.Tribe` object to make an enriched
    **group**-like **output** that contains lists of the :class:`~eqcorrscan.Template`
    objects organized into their groups

    :param tribe: tribe object used to create the inputs to `cluster`
    :type tribe: eqcorrscan.Tribe
    :param groups: output of `cluster` lists of tuples with the
        2nd tuple entry being either the template position in **tribe**
        or the name of the template
        also see :meth:`~eqcorrscan_utils.process.clustering.make_template_list`
    :type groups: list of tuples
    :param savedir: where to save/move outputs from `cluster` and this method, 
        defaults to 'cluster_results'
    :type savedir: str, optional
    :param save_subtribes: should templates be saved to disk in subtribe folders
        within **savedir**? Defaults to False
    :type save_subtribes: bool, optional
    :return:
     - **output** (*dict*) - dictionary with keys corresponding to group numbers
        and values consisting of lists of templates.

        Each of these lists can be directly input to a :class:`~eqcorrscan.Tribe`
        object to reconstitute their functionality via the Tribe class
    """
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    # Save dist_mat.npy
    src = Path().cwd() / 'dist_mat.npy'
    if not os.path.isfile(src):
        Logger.error(f'dist_mat.npy not present in {Path().cwd()}')
    else:
        dest = os.path.join(savedir,'dist_mat.npy')
    os.rename(src, dest)

    # Save cluster_params.csv
    src = Path().cwd() / 'cluster_params.csv'
    if not os.path.isfile(src):
        Logger.error(f'cluster_params.csv not present in {Path().cwd()}')
    else:
        dest = os.path.join(savedir,'cluster_params.csv')
    os.rename(src, dest)

    # replace streams with templates in groups
    output = {}
    for _e, group in enumerate(groups):
        output.update({_e: []})
        for entry in group:
            st, idx = entry
            if isinstance(idx, int):
                output[_e].append(tribe.templates[idx])
            elif isinstance(idx, str):
                output[_e].append(tribe.select(idx))
            else:
                Logger.critical(f'group index {idx} did not look like an index value (int) or name (str)')
        if save_subtribes:
            subtribe_savedir = os.path.join(savedir,f'subtribe_{_e}')
            if not os.path.exists(subtribe_savedir):
                os.makedirs(subtribe_savedir)
            for _t in output[_e]:
                _t.write(os.path.join(subtribe_savedir,f'{_t.name}.tgz'))
    
    # Write grouping to csv
    with open(os.path.join(savedir,'groups.csv'), 'w') as file:
        file.write('template_name,group_id\n')
        for _g, _l in output.items():
            for _t in _l:
                file.write(f'{_t.name},{_g}\n')

    return output




