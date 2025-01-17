import pandas as pd

def reindex_series(ser, ascending=False, inplace=True, start_int=1):
    if not isinstance(ser, pd.Series):
        raise TypeError
    if not isinstance(ascending, bool):
        raise TypeError
    if not isinstance(inplace, bool):
        raise TypeError
    if not isinstance(start_int, int):
        raise TypeError
    # Get a value_counts view of series contents
    _vc = ser.value_counts()
    ngrp = len(_vc)
    if ascending:
        _mapping = {_k: (ngrp - _e) + start_int for _e, _k in enumerate(_vc.index)}
    else:
        _mapping = {_k: _e + start_int for _e, _k in enumerate(_vc.index)}
    if inplace:
        for _k, _v in ser.items():
            if _v in _mapping.keys():
                ser[_k] = _mapping[_v]
        return ser
    else:
        data = []
        for _k, _v in ser.items():
            if _v in _mapping.keys():
                data.append(_mapping[_v])
            else:
                data.append(_v)
        return pd.Series(data=data, index=ser.index.copy(), name=ser.name)
    


def reindex_columns(df, columns, ascending=False, inplace=False, start_int=1):
    if isinstance(df, pd.Series):
        if columns == df.name:
            out = reindex_series(df, ascending=ascending, inplace=inplace, start_int=start_int);
        else:
            raise KeyError
    elif isinstance(df, pd.DataFrame):
        if isinstance(columns, list):
            if all(_e in df.columns for _e in columns):
                pass
            else:
                raise KeyError('Not all specified columns are present in df')
        elif columns in df.columns:
            columns = [columns]
        else:
            raise KeyError(f'column "{columns}" not present in df')
    else:
        raise TypeError
    
    if not inplace:
        out = df.copy()
    else:
        out = df
    for col in columns:
        _ = reindex_series(df[col], ascending=ascending, inplace=True, start_int=start_int);
    
    return out


