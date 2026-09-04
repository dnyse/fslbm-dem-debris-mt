Here's the full 30071 row decoded — every populated field, what it's called, and what it means. Blank fields in the CSV are omitted (they carry no value in the database).

| Group | Characteristic | Value | What it means |
|---|---|---|---|
| **Identification** | ID | 30071 | Unique 5-digit ID; 30000–39999 = Belgium |
| | Institution | ULiège | Team responsible for encoding |
| | Encoder | GM/LB | Initials of the researchers who encoded / checked |
| | Date | 20/10/2022 | Date the encoding was validated |
| **Location** | Type | Railway bridge | Structure type |
| | River | Vesdre | River the bridge crosses |
| | Municipality | Trooz | Municipality it sits in |
| | Structure/street name | Pont ferroviaire Rue Goffontaine | Name of the bridge / street over it |
| | EPSG | 31370 | Coordinate reference system (Belgian Lambert 72) |
| | X reg | 248731 | X coordinate of bridge centre in that CRS [m] |
| | Y reg | 140424 | Y coordinate of bridge centre in that CRS [m] |
| | Lat | 50.56619 | Latitude of bridge centre [°N] |
| | Long | 5.76242 | Longitude of bridge centre [°E] |
| | Curv | 18395 | Curvilinear distance along the river axis from the mouth [m] |
| | Riverbed elevation | 113.0 | Minimum bed elevation below the bridge [m] |
| **Bridge geometry** | Upstream river shape | Straight | Planform of the reach upstream (straight approach flow) |
| | Opening(s) shape | Arched | Shape of the openings water flows through |
| | Width | 10 | Deck dimension *perpendicular to its length* — i.e. streamwise width of the deck [m] |
| | Length | 31.7 | Deck dimension bank to bank — effectively the span across the river [m] |
| | Slope | 0 | Cross-slope of the deck, right→left bank [%] |
| | Angle | 20 | Skew angle between bridge width and the left bank [°] |
| | Thickness | 2.7 | Deck thickness at the opening centre [m] |
| | Elevation | 123.52 | Centre elevation of the deck upper face [m] |
| | River cross-section | Regular | Bed cross-section is (near-)symmetric rectangular/trapezoidal |
| | Abutments | Present | Abutments exist on the banks |
| | Number of piers | 2 | Piers standing in the riverbed (abutments not counted) |
| | Pier width | 2 | Max pier dimension perpendicular to the bed axis [m] |
| | Distance between piers | 8.9-8.9-8.9 | The three opening gaps, left→right bank [m] |
| | Min distance | 8.9 | Smallest gap between adjacent piers/abutments [m] |
| | Max distance | 8.9 | Largest gap between adjacent piers/abutments [m] |
| | Pier shape | Rounded | Shape of the pier nose facing the flow (rounded) |
| | Pier protrusion | 0 | Distance the pier nose projects upstream of the deck [m] — flush here |
| | Handrail material | Metal | Material of the handrail |
| | Handrail porosity | High | Thin elements, large spacing (very open railing) |
| | Structure damage | No | No damage observed after the event |
| **Flood conditions** | Flood event | July 2021 | Name of the flood event |
| | Type of flow | Free surface | Peak water level stayed *below* the deck |
| | Discharge | 535 | Max discharge at the site during the event [m³/s] |
| | Max water elevation | 116.49 | Highest water surface elevation reached [m] |
| | Max water depth | 3.49 | Max elevation minus bed elevation (116.49 − 113.0) [m] |
| **Accumulation** | Clogging | Yes | Debris clogged the opening during the event (certain) |
| | Carpet | Yes | A continuous, compact debris mat formed just upstream |
| | Total length | 10 | Max deposit dimension *parallel* to the river axis [m] |
| | Total width | 12 | Max deposit dimension *normal* to the river axis [m] |
| | Total height | 1 | Max visible vertical dimension of the deposit [m] |
| | Carpet length | 10 | Carpet dimension parallel to river axis [m] |
| | Carpet width | 12 | Carpet dimension normal to river axis [m] |
| | Carpet height | 1 | Visible vertical dimension of the carpet [m] |
| | Volume | 60 | Visible volume of debris blocked at the bridge [m³] |
| | Location at structure | Left bank | Deposit is <80% of span, mainly on the left bank |
| **Main debris content** | Main trunk presence | Yes | A large trunk (vs. opening size) was blocked in the opening |
| | Id main type 1 | 20 → Natural wood | Most abundant element type |
| | Id main type 2 | 27 → Building rubble | Second element type |
| | Id main type 3 | 21 → Anthropogenic wood | Third element type |
| | Volume percentage 1 | 40 | Natural wood share of volume [%] |
| | Volume percentage 2 | 40 | Building rubble share [%] |
| | Volume percentage 3 | 20 | Anthropogenic wood share [%] |

A couple of watch-outs when you pull these into a setup:

The `Width`/`Length` naming is counterintuitive — `Length` (31.7 m) is the cross-stream span and `Width` (10 m) is the streamwise deck depth. For the *deposit*, the axis convention flips: `Total width` (12 m) is cross-stream, `Total length` (10 m) is streamwise (into the flow).

The geometry closes cleanly: 3 × 8.9 + 2 × 2 = 30.7 m against the 31.7 m span, leaving ~1 m for the abutment faces — consistent, unlike 30028. And the debris percentages sum to 100 % (40+40+20), so all of the deposit volume is accounted for by three types with no remainder.