from __future__ import annotations
from dataclasses import dataclass, field
from typing import List
import pandas as pd
from pandas.api.types import is_numeric_dtype
from typing import Dict, Any, Iterator, Tuple, Optional

@dataclass(slots=True)
class FeatureRec:
    id: str                    
    entity: str             
    namespace: str             
    attrs: Dict[str, Any] = field(default_factory=dict)  

    @property
    def properties(self):
        props = self.attrs.copy()
        props["namespace"] = self.namespace
        return props

@dataclass(slots=True)
class ColumnRec:
    id: str                    
    entity: str
    namespace: str
    role: str                  
    attrs: Dict[str, Any] = field(default_factory=dict)

@dataclass(frozen=True, slots=True)
class CrossRef:     
    src: str         
    namespace: str
    target: str

    def to_tuple(self):
        return (self.src,self.namespace,self.target)
    
    def __iter__(self):
        yield self.src
        yield self.namespace
        yield self.target

@dataclass(frozen=True, slots=True)
class OmicData:
    """
    Minimal container for parsed omics data.

    - matrix: numeric measurements
      rows = features; cols = heterogeneous columns (often samples)
    - feature_meta: per-feature metadata (index == matrix.index)
      recommended cols: ["entity","namespace", ...]
    - column_meta: per-column metadata (index == matrix.columns)
      recommended cols: ["entity","namespace","role", ...]
    - cross_ref: long-form cross-refs (axis in {"feature","column"})
      required cols: ["axis","src","namespace","target"]
    """
    name: str
    omics_type: str
    matrix: pd.DataFrame
    feature_meta: List[FeatureRec] = None
    column_meta: List[ColumnRec] = None
    cross_ref: list[CrossRef] | None = None

    def __post_init__(self):
        m = object.__getattribute__(self, "matrix")

        
        if m.index.has_duplicates or m.columns.has_duplicates:
            raise ValueError("Duplicate feature or column IDs")
        if not all(is_numeric_dtype(dt) for dt in m.dtypes):
            raise ValueError("matrix must be numeric-typed")

        fm = object.__getattribute__(self, "feature_meta")
        if fm is not None:
            if not isinstance(fm, list):
                raise ValueError("feature_meta must be a list of FeatureRec")
            ids = [r.id for r in fm]
            if list(m.index) != ids:
                raise ValueError("feature_meta IDs must match matrix.index order")

        
        cm = object.__getattribute__(self, "column_meta")
        if cm is not None:
            if not isinstance(cm, list):
                raise ValueError("column_meta must be a list of ColumnRec")
            ids = [r.id for r in cm]
            if list(m.columns) != ids:
                raise ValueError("column_meta IDs must match matrix.columns order")

    @property
    def features(self) -> pd.Index:
        return self.matrix.index

    @property
    def columns(self) -> pd.Index:
        return self.matrix.columns
    
    def __iter__(self) -> Iterator[Tuple[str, pd.Series, Optional[FeatureRec]]]:
        """
        Yields (feature_id, row_series, feature_meta) for each row in matrix.
        feature_meta is None if feature_meta was not provided.
        """
        fm = self.feature_meta
        for i, fid in enumerate(self.matrix.index):
            meta = fm[i] if fm is not None else None
            yield str(fid), self.matrix.iloc[i, :], meta