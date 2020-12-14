Parameter Options
=================

The parameters to identify the interactions refer to those 
used in `PyPLIF <https://doi.org/10.6026/97320630009325>`_, which 
inspired by the IFP of `Marcou and Rognan <https://doi.org/10.1021/ci600342e>`_. 
These parameters can be modified in ``.hippos/PARAMETERS.py`` in your home
directory after HIPPOS installed.

In general there are two rules that can be modified:

1. Maximum distance value.
2. Interaction angle limit (for hydrogen bond and aromatic interaction).

Here is the content of ``.hippos/PARAMETERS.py`` ::

    '''
        Parameter for interaction distance
        Interaction exist if lower or equal 
        than these values
    '''

    HYDROPHOBIC		= 4.5 
    AROMATIC		= 4.0
    HBOND		= 3.5
    ELECTROSTATIC	= 4.0

    '''
        Parameter for minimum H bond angle.
        Interaction exist if O --- H-D angle
        higher or equal than HBOND_ANGLE value.
    '''

    HBOND_ANGLE     = 135

    '''
        Parameter for aromatic interaction angle.
            Face to Face if:
        AROMATIC_ANGLE_LOW >= Angle between aromatic plane
        Or AROMATIC_ANGLE_HIGH <= Angle between aromatic plane
            Edge to Face if:
        AROMATIC_ANGLE_LOW < Angle between aromatic plane < 
        AROMATIC_ANGLE_HIGH
    '''

    AROMATIC_ANGLE_LOW  = 30.0
    AROMATIC_ANGLE_HIGH = 150.0


The explanation of each value can be seen from the comment above the
value assignment.

.. _simplified-rule:

Simplified Bitstring Rule
-------------------------

HIPPOS provides simplified interaction bitstring to reduce unnecessary
string in the output. The bitstring simplification rule can be seen in
the table below. Number 1-7 in the header of the table correspond to
the interaction type described in PyPLIF

.. code-block::

    1 Correspond to Hydrophobic
    2 Correspond to Aromatic Face to Face
    3 Correspond to Aromatic Edge to Face
    4 Correspond to H-bond (protein is donor)
    5 Correspond to H-bond (protein is acceptor)
    6 Correspond to Electrostatic (protein +)
    7 Correspond to Electrostatic (protein -)

.. raw:: html
    
    <embed>
    <style type="text/css">
    .tg  {border:none;border-collapse:collapse;border-color:#9ABAD9;border-spacing:0;}
    .tg td{background-color:#EBF5FF;border-color:#9ABAD9;border-style:solid;border-width:0px;color:#444;
      font-family:Arial, sans-serif;font-size:14px;overflow:hidden;padding:5px 15px;word-break:normal;}
    .tg th{background-color:#409cff;border-color:#9ABAD9;border-style:solid;border-width:0px;color:#fff;
      font-family:Arial, sans-serif;font-size:14px;font-weight:bold;overflow:hidden;padding:5px 15px;word-break:normal;}
    .tg .tg-h1fx{background-color:#D2E4FC;border-color:inherit;font-family:"Courier New", Courier, monospace !important;;
      text-align:center;vertical-align:top}
    .tg .tg-juju{font-family:"Courier New", Courier, monospace !important;;text-align:left;vertical-align:top}
    .tg .tg-3ib7{border-color:inherit;font-family:"Courier New", Courier, monospace !important;;text-align:center;vertical-align:top}
    .tg .tg-64ye{background-color:#D2E4FC;font-family:"Courier New", Courier, monospace !important;;text-align:left;vertical-align:top}
    </style>
    <table class="tg">
    <thead>
      <tr>
        <th class="tg-3ib7">AA</th>
        <th class="tg-3ib7">1</th>
        <th class="tg-3ib7">2</th>
        <th class="tg-juju">3</th>
        <th class="tg-3ib7">4</th>
        <th class="tg-3ib7">5</th>
        <th class="tg-3ib7">6</th>
        <th class="tg-3ib7">7</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td class="tg-h1fx">ALA</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">CYS</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;<br></td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">ASP</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
      </tr>
      <tr>
        <td class="tg-3ib7">GLU</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">&#10003;</td>
      </tr>
      <tr>
        <td class="tg-h1fx">PHE</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-64ye">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">GLY</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">HIS</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-64ye">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">ILE</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">LYS</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">LEU</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">MET</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">ASN</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;<br></td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">PRO</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">GLN</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">ARG</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;<br></td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">SER</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-<br></td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">THR</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">VAL</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
      <tr>
        <td class="tg-h1fx">TRP</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-64ye">-</td>
        <td class="tg-h1fx">&#10003;</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
        <td class="tg-h1fx">-</td>
      </tr>
      <tr>
        <td class="tg-3ib7">TYR</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-juju">-</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">&#10003;</td>
        <td class="tg-3ib7">-</td>
        <td class="tg-3ib7">-</td>
      </tr>
    </tbody>
    </table>
    </br>
    </embed>